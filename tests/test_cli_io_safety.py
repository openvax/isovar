"""Every command checks its outputs before working and never overwrites its inputs.

Each case comes from running the CLI on Sid data (#661): outputs replaced the
VCF or input JSON, failed after the whole analysis, or ignored ``.gz``.
"""

from argparse import ArgumentParser, Namespace
import gzip
import json

import pytest

from isovar.cli.commands import run
from isovar.cli.validation import CommandInputError, check_outputs, input_error
from isovar.json_files import json_stem, read_json, write_json


def test_an_output_may_not_overwrite_an_input_named_in_any_option(tmp_path):
    vcf, maf, profile = (tmp_path / name for name in ("calls.vcf", "calls.maf", "profile.json"))
    for path in (vcf, maf, profile):
        path.write_text("keep")
    args = Namespace(vcf=[str(vcf)], maf=str(maf), read_end_profile=profile, output=None, count=3)
    for target in (vcf, maf, profile):
        args.output = str(target)
        with pytest.raises(CommandInputError, match="would overwrite an input file"):
            check_outputs(args, "output")
    args.output = str(tmp_path / "new.json")
    check_outputs(args, "output")
    # A previous output may be replaced.
    (tmp_path / "new.json").write_text("old")
    check_outputs(args, "output")
    assert all(path.read_text() == "keep" for path in (vcf, maf, profile))


def test_outputs_are_files_in_existing_directories_and_distinct(tmp_path):
    (tmp_path / "a-directory").mkdir()
    (tmp_path / "a-file").write_text("")
    args = Namespace(output=str(tmp_path / "a-directory"), plot_dir=str(tmp_path / "a-file"))
    with pytest.raises(CommandInputError, match="Output is a directory"):
        check_outputs(args, "output")
    with pytest.raises(CommandInputError, match="Output directory is an existing file"):
        check_outputs(args, directories=("plot_dir",))
    with pytest.raises(CommandInputError, match="Output directory does not exist"):
        check_outputs(args, str(tmp_path / "missing" / "out.json"))
    with pytest.raises(CommandInputError, match="must be distinct"):
        check_outputs(args, str(tmp_path / "x.json"), str(tmp_path / "x.json"))


def test_input_errors_are_one_line_without_the_usage_block(capsys):
    parser = ArgumentParser(prog="isovar example")
    parser.add_argument("--many", help="an option whose usage would bury the message")
    with pytest.raises(SystemExit) as error:
        input_error(parser, CommandInputError("Output x.json would overwrite an input file"))
    assert error.value.code == 2
    assert capsys.readouterr().err == "isovar example: error: Output x.json would overwrite an input file\n"


def test_json_files_compress_by_name_and_replace_atomically(tmp_path, monkeypatch):
    for name in ("result.json", "result.json.gz"):
        path = tmp_path / name
        write_json(path, {"b": [1, 2]})
        assert read_json(path) == {"b": [1, 2]}
        assert (path.read_bytes()[:2] == b"\x1f\x8b") == name.endswith(".gz")
    path = tmp_path / "result.json"

    def fail(*args, **kwargs):
        raise TypeError("not serializable")
    monkeypatch.setattr(json, "dumps", fail)
    with pytest.raises(TypeError):
        write_json(path, {"replaced": True})
    assert read_json(path) == {"b": [1, 2]}
    assert sorted(p.name for p in tmp_path.iterdir()) == ["result.json", "result.json.gz"]
    assert json_stem(tmp_path / "x.json.gz").name == json_stem(tmp_path / "x.json").name == "x"


def test_protein_hypotheses_writes_gzip_with_the_tsv_beside_the_stem(tmp_path):
    from isovar.protein_hypotheses import protein_hypotheses_tsv_path, write_protein_hypotheses

    assert protein_hypotheses_tsv_path(tmp_path / "x.json.gz") == tmp_path / "x.tsv"
    assert protein_hypotheses_tsv_path(tmp_path / "x.json") == tmp_path / "x.tsv"
    assert protein_hypotheses_tsv_path(tmp_path / "x.out") == tmp_path / "x.out.tsv"
    export = dict(schema="isovar.protein_hypotheses.v2", events=[])
    paths = write_protein_hypotheses(export, tmp_path / "x.json.gz")
    with gzip.open(paths["json"], "rt") as stream:
        assert json.load(stream)["schema"] == "isovar.protein_hypotheses.v2"
    assert b"\r\n" not in paths["tsv"].read_bytes()


@pytest.mark.parametrize("command, outputs", [
    ("protein-hypotheses", ["--output", "{vcf}"]),
    ("patient-background-pairs", ["--output", "{vcf}"]),
    ("run", ["--output", "{directory}"]),
])
def test_commands_refuse_outputs_before_reading_reads(command, outputs, tmp_path, capsys):
    vcf = tmp_path / "calls.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    (tmp_path / "directory").mkdir()
    values = [value.format(vcf=vcf, directory=tmp_path / "directory") for value in outputs]
    with pytest.raises(SystemExit) as error:
        run([command, "--vcf", str(vcf), "--bam", str(tmp_path / "not-read.bam"), "--genome", "GRCh38",
             "--sample-id", "s"][:None if command != "run" else -2] + values)
    assert error.value.code == 2
    err = capsys.readouterr().err
    assert err.startswith("isovar %s: error: Output" % command) and "usage:" not in err
    assert vcf.read_text().startswith("##fileformat")


def test_json_commands_read_and_write_gzip(tmp_path, capsys):
    interpretations = tmp_path / "in.json.gz"
    write_json(interpretations, [{"not": "an object"}])
    with pytest.raises(SystemExit) as error:
        run(["allele-interpretations", "--input", str(interpretations), "--bam", str(tmp_path / "x.bam"),
             "--sample-id", "s", "--output", str(interpretations)])
    assert error.value.code == 2 and "would overwrite an input file" in capsys.readouterr().err
    with pytest.raises(SystemExit) as error:
        run(["allele-interpretations", "--input", str(interpretations), "--bam", str(tmp_path / "x.bam"),
             "--sample-id", "s", "--output", str(tmp_path / "out.json.gz")])
    assert error.value.code == 2 and "exactly the keys" in capsys.readouterr().err


def test_antigen_report_counts_blank_calls_as_unmeasured_and_guards_its_inputs(tmp_path, capsys):
    from isovar.antigen_report import load_small_variants, small_variant_summary

    tsv = tmp_path / "calls.tsv"
    fields = ["sample_id", "call_id", "ref", "alt", "exon_overlap", "alt_reads", "protein_count",
              "orfdata_coding_block_overlap"]
    rows = [["T1", "a", "G", "T", "True", "7", "1", "False"],
            # An unresolved chromosome: the evidence report leaves counts blank.
            ["T1", "b", "C", "A", "True", "", "", "False"]]
    tsv.write_text("\n".join("\t".join(row) for row in [fields] + rows) + "\n")
    summary = small_variant_summary(load_small_variants(tsv, "T1"))
    assert summary["rna_unmeasured"] == 1 and summary["rna_alt_observed"] == 1
    assert summary["exonic_reads_by_class"] == {"SNV": [7]}
    with pytest.raises(SystemExit) as error:
        run(["antigen-report", "--sample-id", "T1", "--output", str(tmp_path / "r.pdf"),
             "--small-variant-tsv", str(tsv), "--summary-json", str(tsv)])
    assert error.value.code == 2 and "would overwrite an input file" in capsys.readouterr().err
    assert tsv.read_text().startswith("sample_id")


def test_cohort_reports_are_never_rendered_from_a_source_selection(tmp_path):
    from isovar.cohort import run_cohort

    with pytest.raises(ValueError, match="Reports cover every source"):
        run_cohort(tmp_path / "manifest.json", tmp_path / "out", source_ids=["T1"])
    with pytest.raises(SystemExit) as error:
        run(["cohort", "--manifest", str(tmp_path / "manifest.json"), "--output-dir", str(tmp_path / "out"),
             "--render-only", "--resume"])
    assert error.value.code == 2


def test_sv_source_rejects_an_unknown_assembly_before_reading_rows(tmp_path, capsys):
    from isovar.sv_rna_audit import file_sha256

    vcf = tmp_path / "sv.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    with pytest.raises(SystemExit) as error:
        run(["sv-source", "--vcf", str(vcf), "--sha256", file_sha256(vcf), "--reference-name", "banana",
             "--output", str(tmp_path / "out.json")])
    assert error.value.code == 2
    assert "Unknown reference assembly 'banana'" in capsys.readouterr().err
    assert not (tmp_path / "out.json").exists()


def test_a_failed_isoform_reconstruction_leaves_no_output_directory(tmp_path, capsys):
    # The directory used to be created first, so repeating a failed run hit "File exists".
    output = tmp_path / "isoforms"
    command = ["isoforms", "reconstruct", "--bam", str(tmp_path / "missing.bam"), "--sample-id", "s",
               "--gene", "TPST1", "--genome-fasta", str(tmp_path / "missing.fa"),
               "--gtf", str(tmp_path / "missing.gtf"), "--output-dir", str(output)]
    for _ in range(2):
        with pytest.raises(SystemExit) as error:
            run(command)
        assert error.value.code == 2 and "File exists" not in capsys.readouterr().err
        assert not output.exists()
