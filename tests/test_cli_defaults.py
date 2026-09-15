"""CLI adapters must preserve Python defaults and explicit user choices."""

from inspect import signature

import pytest

from isovar import ProteinSequenceCreator, ReadCollector, run_isovar
from isovar import default_parameters as defaults
from isovar.cli import rna_args
from isovar.cli.filter_args import filter_threshold_dict_from_args
from isovar.cli.main_args import make_isovar_arg_parser
from isovar.cli.protein_sequence_args import (
    make_protein_sequences_arg_parser,
    protein_sequence_creator_from_args,
)
from isovar.cli.reference_context_args import make_reference_context_arg_parser
from isovar.cli.translation_args import (
    make_translation_arg_parser,
    protein_sequence_creator_kwargs_from_args,
)
from isovar.cli.variant_sequences_args import make_variant_sequences_arg_parser
from isovar.reference_context_helpers import reference_contexts_for_variant, reference_contexts_generator
from isovar.variant_sequence_creator import VariantSequenceCreator

from .mock_objects import MockAlignmentFile, make_pysam_read


@pytest.mark.parametrize("factory", [
    make_isovar_arg_parser, make_protein_sequences_arg_parser, make_translation_arg_parser,
])
def test_cli_protein_defaults_match_python(factory):
    args = factory().parse_args(["--bam", "unused.bam"])
    cli = ProteinSequenceCreator(**protein_sequence_creator_kwargs_from_args(args))
    api = ProteinSequenceCreator()
    for name in signature(ProteinSequenceCreator).parameters:
        assert getattr(cli, name) == getattr(api, name), name
    assert vars(cli._variant_sequence_creator) == vars(api._variant_sequence_creator)


def test_cli_rna_defaults_match_python_and_merge_mates():
    args = make_isovar_arg_parser().parse_args(["--bam", "unused.bam"])
    cli = rna_args.read_collector_from_args(args)
    api = ReadCollector()
    assert vars(cli) == vars(api)
    reads = [
        make_pysam_read("ACCGTG", "6M", name="pair", reference_start=0),
        make_pysam_read("CGTGAA", "6M", name="pair", reference_start=2),
    ]
    for collector in (cli, api):
        merged = collector.get_locus_reads(MockAlignmentFile(["1"], reads), "1", 3, 4)
        assert len(merged) == 1
        assert merged[0].sequence == "ACCGTGAA"
        assert merged[0].source_read_count == 2


def test_explicit_rna_options_and_legacy_namespace():
    args = make_isovar_arg_parser().parse_args([
        "--bam", "unused.bam", "--min-mapping-quality", "17",
        "--use-duplicate-reads", "--drop-secondary-alignments",
        "--use-soft-clipped-bases", "--no-merge-overlapping-fragments",
    ])
    assert vars(rna_args.read_collector_from_args(args)) == vars(ReadCollector(
        min_mapping_quality=17, use_duplicate_reads=True, use_secondary_alignments=False,
        use_soft_clipped_bases=True, merge_overlapping_fragments=False))
    del args.merge_overlapping_fragments
    assert rna_args.read_collector_from_args(args).merge_overlapping_fragments == \
        defaults.MERGE_OVERLAPPING_FRAGMENTS


def test_cli_variant_sequence_defaults_match_python():
    args = make_variant_sequences_arg_parser(add_sequence_length_arg=True).parse_args(["--bam", "unused"])
    creator = VariantSequenceCreator()
    assert args.variant_sequence_length == creator.preferred_sequence_length
    assert args.min_variant_sequence_coverage == creator.min_variant_sequence_coverage
    assert args.variant_sequence_assembly == creator.variant_sequence_assembly


def test_filter_defaults_match_python_and_overrides_do_not_remove_other_filters():
    parser = make_isovar_arg_parser()
    args = parser.parse_args(["--bam", "unused"])
    assert filter_threshold_dict_from_args(args) == defaults.DEFAULT_FILTER_THRESHOLDS
    args = parser.parse_args([
        "--bam", "unused", "--min-alt-rna-reads", "7", "--min-alt-rna-fragments", "4",
        "--min-alt-rna-fraction", "0.1", "--min-ratio-alt-to-other-fragments", "5",
    ])
    expected = defaults.DEFAULT_FILTER_THRESHOLDS.copy()
    expected.update(min_num_alt_reads=7, min_num_alt_fragments=4,
                    min_fraction_alt_fragments=0.1, min_ratio_alt_to_other_fragments=5)
    assert filter_threshold_dict_from_args(args) == expected
    expected["max_num_other_reads"] = -1
    assert defaults.DEFAULT_FILTER_THRESHOLDS["max_num_other_reads"] != -1


def test_main_keeps_public_default_filter_imports():
    from isovar.main import DEFAULT_FILTER_THRESHOLDS, DEFAULT_FILTER_FLAGS
    assert DEFAULT_FILTER_THRESHOLDS is defaults.DEFAULT_FILTER_THRESHOLDS
    assert DEFAULT_FILTER_FLAGS is defaults.DEFAULT_FILTER_FLAGS


@pytest.mark.parametrize("threads", [None, 4])
def test_cli_and_api_pass_same_decompression_threads(monkeypatch, threads):
    import isovar.main as main
    calls = []

    def alignment_file(path, **kwargs):
        calls.append((path, kwargs))
        return MockAlignmentFile([], [])

    monkeypatch.setattr(main, "AlignmentFile", alignment_file)
    monkeypatch.setattr(rna_args, "AlignmentFile", alignment_file)
    cli_options = ["--bam", "unused.bam"]
    api_options = {}
    if threads is not None:
        cli_options += ["--num-rna-decompression-threads", str(threads)]
        api_options["decompression_threads"] = threads
    run_isovar([], "unused.bam", **api_options)
    rna_args.alignment_file_from_args(make_isovar_arg_parser().parse_args(cli_options))
    assert calls[0] == calls[1]
    assert calls[0][1]["threads"] == (defaults.NUM_RNA_DECOMPRESSION_THREADS if threads is None else threads)


def test_explicit_protein_options_and_vaxrank_peptide_override():
    args = make_protein_sequences_arg_parser().parse_args([
        "--bam", "unused", "--protein-sequence-length", "60",
        "--protein-context-peptide-length", "15", "--protein-sequence-preference", "support",
        "--min-variant-sequence-coverage", "3", "--min-transcript-prefix-length", "12",
        "--max-reference-transcript-mismatches", "4", "--count-mismatches-after-variant",
        "--disable-variant-sequence-assembly", "--max-protein-sequences-per-variant", "0",
        "--min-protein-sequence-support-fraction", "0.9",
    ])
    cli = protein_sequence_creator_from_args(args)
    api = ProteinSequenceCreator(
        protein_sequence_length=60, protein_context_peptide_length=15,
        protein_sequence_preference="support", min_variant_sequence_coverage=3,
        min_transcript_prefix_length=12, max_transcript_mismatches=4,
        count_mismatches_after_variant=True, variant_sequence_assembly=False,
        max_protein_sequences_per_variant=0, min_protein_sequence_support_fraction=0.9)
    for name in signature(ProteinSequenceCreator).parameters:
        assert getattr(cli, name) == getattr(api, name), name

    args.protein_sequence_length = None
    args.protein_context_peptide_length = None
    args.vaccine_peptide_length = 31
    cli = protein_sequence_creator_from_args(args)
    assert cli.protein_context_peptide_length == 31
    assert cli.protein_sequence_length == 61
    args.protein_context_peptide_length = 15
    assert protein_sequence_creator_from_args(args).protein_sequence_length == 29


@pytest.mark.parametrize("factory", [
    make_isovar_arg_parser, make_protein_sequences_arg_parser, make_translation_arg_parser,
])
def test_protein_commands_reject_ignored_reference_context_option(factory, capsys):
    parser = factory()
    assert "--reference-context-size" not in parser.format_help()
    with pytest.raises(SystemExit) as error:
        parser.parse_args(["--bam", "unused", "--reference-context-size", "5"])
    assert error.value.code == 2
    assert "unrecognized arguments: --reference-context-size" in capsys.readouterr().err


@pytest.mark.parametrize("value", ["0", "-7"])
def test_reference_context_command_rejects_nonpositive_sizes(value, capsys):
    with pytest.raises(SystemExit) as error:
        make_reference_context_arg_parser().parse_args(["--reference-context-size", value])
    assert error.value.code == 2
    assert "positive integer" in capsys.readouterr().err


def test_reference_context_command_passes_size_to_api(monkeypatch):
    from isovar.cli import reference_context_args
    parser = make_reference_context_arg_parser()
    default = parser.parse_args([]).reference_context_size
    assert default == defaults.REFERENCE_CONTEXT_SIZE
    for fn in (reference_contexts_for_variant, reference_contexts_generator):
        assert signature(fn).parameters["context_size"].default == default
    monkeypatch.setattr(reference_context_args, "variant_collection_from_args", lambda args: ["variant"])
    calls = []

    def contexts(**kwargs):
        calls.append(kwargs)
        return []

    monkeypatch.setattr(reference_context_args, "reference_contexts_generator", contexts)
    monkeypatch.setattr(reference_context_args, "variants_to_reference_contexts_dataframe", list)
    reference_context_args.reference_contexts_dataframe_from_args(
        parser.parse_args(["--reference-context-size", "7"]))
    assert calls == [{"variants": ["variant"], "context_size": 7}]
