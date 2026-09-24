"""Explicit fusion sequence/frame validation, independent of small variants."""
from dataclasses import asdict, replace
from inspect import signature
import gzip
import hashlib
import json
from pathlib import Path

import pytest
import pysam

from isovar import FusionBlock, FusionBreakpoint, FusionRead, FusionReference, FusionTranscript, reconstruct_fusion
from isovar.cli import commands
from isovar.cli.isovar_fusion import make_parser
from isovar.fusion import fusion_from_dict
from tests.testing_helpers import fusion_input
from .data.osteosarc.expansion.references import translate


def example(donor_strand="+", acceptor_strand="+", cut=12, insert="", partial=0):
    donor_sequence = "ATGGCTGCTAAACCTGAATAA"
    acceptor_sequence = "ATGGGTCCCTTTGAATAA"
    donor = FusionReference("donor.1", "test", "test annotation 1", "d", donor_strand,
                            ((100, 100 + len(donor_sequence)),), donor_sequence, 0, len(donor_sequence))
    acceptor = FusionReference("acceptor.1", "test", "test annotation 1", "a", acceptor_strand,
                               ((200, 200 + len(acceptor_sequence)),), acceptor_sequence, 0, len(acceptor_sequence))
    sequence = donor_sequence[partial:cut] + insert + acceptor_sequence[3:]
    j0, j1 = cut - partial, cut - partial + len(insert)
    de = 100 + cut if donor_strand == "+" else 100 + len(donor_sequence) - cut
    ae = 203 if acceptor_strand == "+" else 200 + len(acceptor_sequence) - 3
    blocks = (FusionBlock(0, j0, "d", 100 + partial if donor_strand == "+" else de,
                          de if donor_strand == "+" else 100 + len(donor_sequence) - partial, donor_strand),
              FusionBlock(j1, len(sequence), "a", ae if acceptor_strand == "+" else 200,
                          200 + len(acceptor_sequence) if acceptor_strand == "+" else ae, acceptor_strand))
    fusion = FusionTranscript("donor--acceptor", "test", sequence, j0, j1,
        FusionBreakpoint("d", de, donor_strand), FusionBreakpoint("a", ae, acceptor_strand), blocks,
        dict(sample_id="sample", method="test assembly", version="1", parameters={}, source="synthetic", contig_id="contig"))
    reads = tuple(FusionRead("sample", "library", "fragment" + str(i), "read" + str(i),
                            "synthetic", 0, 0, sequence, blocks) for i in range(2))
    return fusion, (donor, acceptor), reads


@pytest.mark.parametrize("donor_strand", ["+", "-"])
@pytest.mark.parametrize("acceptor_strand", ["+", "-"])
@pytest.mark.parametrize("cut,insert", [(12, ""), (13, ""), (12, "A"), (13, "TG"), (12, "GCT")])
def test_orientations_split_codons_insertions_and_independent_translation(donor_strand, acceptor_strand, cut, insert):
    fusion, refs, reads = example(donor_strand, acceptor_strand, cut, insert)
    result = reconstruct_fusion(fusion, refs, reads, peptide_lengths=(2, 3))
    assert result["status"] == "translated"
    protein = result["paths"][0]["translations"][0]
    expected, stop = translate(fusion.sequence, annotated_start=True)
    assert protein["amino_acids"] == expected
    assert protein["ends_with_stop_codon"] == stop
    assert protein["junction_in_translated_cds"] == [cut, cut + len(insert)]
    assert protein["acceptor_frames"][0]["in_frame"] == ((cut + len(insert)) % 3 == 0)
    for peptide in protein["candidate_peptides"]:
        lo, hi = peptide["protein_interval"]
        assert peptide["sequence"] == translate(fusion.sequence[3 * lo:3 * hi], annotated_start=lo == 0)[0]
        assert any(3 * lo < b < 3 * hi for b in (cut, cut + len(insert)))
    for peptide in protein["downstream_frameshift_peptides"]:
        lo, hi = peptide["protein_interval"]
        assert 3 * lo >= cut + len(insert)
        assert peptide["sequence"] == expected[lo:hi]
        assert not protein["acceptor_frames"][0]["in_frame"]


@pytest.mark.parametrize("partial", [1, 2, 3, 4, 5, 6])
def test_partial_cds_is_conditional_not_invented_start(partial):
    fusion, refs, reads = example(partial=partial)
    result = reconstruct_fusion(fusion, refs, reads, peptide_lengths=(2,))
    protein = result["paths"][0]["translations"][0]
    start = (-partial) % 3
    assert not protein["complete_5prime"] and protein["cds_start"] is None
    assert protein["translation_start"] == start
    assert protein["frame_evidence"][0]["upstream_frame_assumed"]
    assert protein["amino_acids"] == translate(fusion.sequence[start:])[0]
    assert protein["junction_in_translated_cds"] == [12 - partial - start] * 2


def test_compatible_transcripts_are_a_set_not_a_best_isoform():
    fusion, refs, reads = example()
    refs += (replace(refs[0], transcript_id="donor.2"),)
    result = reconstruct_fusion(fusion, refs, reads)
    assert result["status"] == "translated"
    assert result["paths"][0]["translations"][0]["transcript_ids"] == ["donor.1", "donor.2"]
    result = reconstruct_fusion(fusion, refs + (replace(refs[0], transcript_id="noncoding", cds_start=None, cds_end=None),), reads)
    assert result["status"] == "ambiguous"
    assert result["reasons"] == ["donor_CDS_unavailable:noncoding"]


def test_alternative_upstream_frames_are_not_ranked():
    fusion, refs, reads = example(partial=3, cut=9)
    # Identical observed donor interval, different annotated upstream exon/start.
    alternative = FusionReference("donor.other_frame", "test", "test annotation 1", "d", "+",
        ((0, 4), (103, 109), (300, 305)), "ATGA" + refs[0].sequence[3:9] + "CCTAA", 0, 15)
    result = reconstruct_fusion(fusion, refs + (alternative,), reads)
    assert result["status"] == "ambiguous"
    assert {p["translation_start"] for p in result["paths"][0]["translations"]} == {0, 2}
    assert all(p["cds_start"] is None for p in result["paths"][0]["translations"])


def test_no_cds_no_reference_or_inconsistent_donor_never_uses_longest_orf():
    fusion, refs, reads = example()
    assert reconstruct_fusion(fusion, (), reads)["status"] == "unresolved"
    noncoding = replace(refs[0], cds_start=None, cds_end=None)
    assert reconstruct_fusion(fusion, (noncoding,), reads)["paths"][0]["translations"] == []
    changed = replace(fusion, sequence=fusion.sequence[:6] + "A" + fusion.sequence[7:])
    changed_reads = tuple(replace(r, sequence=changed.sequence) for r in reads)
    assert reconstruct_fusion(changed, refs, changed_reads)["status"] == "unresolved"


def test_utr_join_early_stop_and_no_acceptor_codon_are_explicit():
    fusion, refs, reads = example()
    utr = replace(refs[0], sequence=refs[0].sequence[:12] + "ATGAAATAA", cds_start=12, cds_end=21)
    result = reconstruct_fusion(fusion, (utr,), reads)
    assert result["reasons"] == ["junction_before_donor_CDS:donor.1"]
    fusion, refs, reads = example(insert="TAA")
    result = reconstruct_fusion(fusion, refs, reads)
    assert result["paths"][0]["translations"] == []
    assert result["reasons"] == ["no_translated_acceptor_context:donor.1"]
    fusion, refs, reads = example(cut=20)
    assert reconstruct_fusion(fusion, refs, reads)["reasons"] == ["junction_after_donor_CDS:donor.1"]


def test_incomplete_donor_codon_and_disrupted_start_do_not_claim_a_fusion_protein():
    for partial in (10, 11):
        result = reconstruct_fusion(*example(partial=partial))
        assert result["paths"][0]["translations"] == []
        assert result["reasons"] == ["no_translated_donor_context:donor.1"]
    result = reconstruct_fusion(*example(cut=1))
    assert result["paths"][0]["translations"] == []
    assert result["reasons"] == ["fusion_disrupts_annotated_start:donor.1"]


def test_duplicate_products_supplementary_records_and_mates_do_not_inflate_fragments():
    fusion, refs, reads = example()
    result = reconstruct_fusion(fusion, refs, reads + (replace(reads[0], source="processed-copy"),))
    assert result["evidence"]["reads"] == result["evidence"]["direct_fragments"] == 2
    mates = (reads[0], replace(reads[1], fragment_id=reads[0].fragment_id))
    result = reconstruct_fusion(fusion, refs, mates)
    assert result["status"] == "insufficient_support"
    assert result["evidence"]["reads"] == 2 and result["evidence"]["fragments"] == 1
    with pytest.raises(ValueError, match="Conflicting observations"):
        reconstruct_fusion(fusion, refs, reads + (replace(reads[0], source_query_start=1),))


def test_paired_mates_share_query_name_without_becoming_duplicates():
    fusion, refs, reads = example()
    first = replace(reads[0], read_id="paired", fragment_id="template", mate_number=1)
    second = replace(reads[1], read_id="paired", fragment_id="template", mate_number=2)

    result = reconstruct_fusion(fusion, refs, (first, second))

    assert result["evidence"]["reads"] == 2
    assert result["evidence"]["fragments"] == 1
    assert result["evidence"]["direct_fragments"] == 1
    assert result["status"] == "insufficient_support"
    duplicate_first = replace(first, source="processed-copy")
    deduplicated = reconstruct_fusion(
        fusion, refs, (first, duplicate_first, second)
    )
    assert deduplicated["evidence"]["reads"] == 2
    with pytest.raises(ValueError, match="Conflicting observations"):
        reconstruct_fusion(
            fusion,
            refs,
            (first, replace(first, source_query_start=1), second),
        )


def test_unmapped_or_unobserved_sequence_does_not_establish_support():
    fusion, refs, reads = example()
    assert reconstruct_fusion(fusion, refs)["status"] == "insufficient_support"
    assert reconstruct_fusion(fusion, refs, tuple(replace(r, blocks=()) for r in reads))["status"] == "insufficient_support"
    # A junction witness alone does not validate unobserved contig flanks.
    short = tuple(replace(r, sequence=r.sequence[:-3], blocks=(r.blocks[0],
        replace(r.blocks[1], query_end=r.blocks[1].query_end-3, reference_end=r.blocks[1].reference_end-3))) for r in reads)
    assert reconstruct_fusion(fusion, refs, short)["reasons"] == ["unobserved_fusion_sequence"]


@pytest.mark.parametrize("field,value,error", [
    ("sample_id", "other", "another sample"), ("sequence", "AAAA", "does not match"),
    ("cdna_start", 1, "does not match"),
])
def test_conflicting_read_evidence_is_rejected(field, value, error):
    fusion, refs, reads = example()
    with pytest.raises(ValueError, match=error):
        reconstruct_fusion(fusion, refs, (replace(reads[0], **{field:value}),))


def test_mapping_validation():
    fusion, refs, reads = example()
    with pytest.raises(ValueError, match="breakpoints"):
        replace(fusion, donor=replace(fusion.donor, position=fusion.donor.position+1))
    with pytest.raises(ValueError, match="orientation"):
        replace(fusion, donor=replace(fusion.donor, strand="-"))
    with pytest.raises(ValueError, match="conflicting partner"):
        reconstruct_fusion(fusion, refs, (replace(reads[0], blocks=(replace(reads[0].blocks[0],contig="other"),)),))
    with pytest.raises(ValueError, match="Overlapping or alternative"):
        reconstruct_fusion(fusion, refs, (replace(reads[0], blocks=reads[0].blocks*2),))
    with pytest.raises(ValueError, match="overlap"):
        replace(fusion, blocks=fusion.blocks*2)
    with pytest.raises(ValueError, match="Duplicate reference"):
        reconstruct_fusion(fusion, refs*2, reads)


@pytest.mark.parametrize("change", [dict(sequence="N"), dict(junction_start=True), dict(junction_end=999),
                                     dict(provenance={}), dict(blocks=())])
def test_invalid_fusion_inputs(change):
    with pytest.raises(ValueError):
        replace(example()[0], **change)


@pytest.mark.parametrize("change", [dict(cds_start=1), dict(cds_end=20), dict(genetic_code=2),
                                     dict(exons=((100,110),)), dict(strand="?")])
def test_invalid_reference_inputs(change):
    with pytest.raises(ValueError):
        replace(example()[1][0], **change)


@pytest.mark.parametrize("mate_number", [True, -1, 3])
def test_invalid_mate_number(mate_number):
    with pytest.raises(ValueError, match="mate_number"):
        replace(example()[2][0], mate_number=mate_number)


@pytest.mark.parametrize("options", [dict(min_fragments=0), dict(min_fragments=True),
                                      dict(peptide_lengths=[]), dict(peptide_lengths=[0])])
def test_invalid_parameters(options):
    with pytest.raises(ValueError):
        reconstruct_fusion(*example(), **options)


def test_json_cli_defaults_and_round_trip(tmp_path):
    fusion, refs, reads = example()
    data = dict(fusion=asdict(fusion), references=[asdict(r) for r in refs], reads=[asdict(r) for r in reads])
    source, output = tmp_path / "input.json", tmp_path / "output.json"
    source.write_text(json.dumps(data))
    assert fusion_from_dict(json.loads(source.read_text())) == (fusion, refs, reads)
    args = make_parser().parse_args(["--input", str(source), "--output", str(output)])
    for key in ("peptide_lengths", "min_fragments"):
        assert getattr(args,key) == signature(reconstruct_fusion).parameters[key].default
    commands.run(["fusion", "--input", str(source), "--output", str(output)])
    assert json.loads(output.read_text()) == json.loads(json.dumps(reconstruct_fusion(fusion, refs, reads)))
    source.write_text('{}')
    with pytest.raises(SystemExit) as error:
        commands.run(["fusion", "--input", str(source), "--output", str(output)])
    assert error.value.code == 2


def test_original_osteosarc_rna_windows_remain_unresolved_not_fake_proteins():
    corpus = Path(__file__).parent / "data/fusions/corpus"
    examples = json.loads((corpus / "manifest.json").read_text())
    checked = []
    for entry in examples:
        if "input" not in entry:
            assert entry["event"].startswith("PARD3B") and entry["sample"] == "T2"
            continue
        path = corpus / entry["input"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == entry["sha256"]
        data = json.loads(gzip.decompress(path.read_bytes()))
        fusion, refs, reads = fusion_from_dict(fusion_input(data))
        result = reconstruct_fusion(fusion, refs, reads)
        assert result["status"] == "unresolved"
        assert result["paths"][0]["translations"] == []
        assert result["evidence"]["direct_fragments"] == entry["fragments"]
        if fusion.event_id.startswith("TPST1"):
            assert "junction_before_donor_CDS:ENST00000304842" in result["reasons"]
        else:
            assert result["reasons"] == ["no_exact_collinear_annotated_donor"]
        if fusion.event_id.startswith("PARD3B"):
            assert fusion.sequence[fusion.junction_start:fusion.junction_end] == "TGTGCTGATATG"
        header = pysam.AlignmentHeader.from_references(["chr" + str(i) for i in range(1, 23)], [300000000] * 22)
        for observation, original in zip(reads, data["original_records"]):
            record = pysam.AlignedSegment.fromstring(original["sam"], header)
            partner = pysam.AlignedSegment.fromstring(original["partner_sam"], header)
            from tests.data.fusions.build_osteosarc import extract
            reconstructed = extract(record, fusion.donor.contig, fusion.donor.position,
                                    fusion.acceptor.contig, fusion.acceptor.position+1,
                                    fusion.junction_start, [record, partner])
            assert reconstructed is not None
            assert reconstructed['sequence'] == fusion.sequence
            assert reconstructed['blocks'] == data['fusion']['blocks']
            lo = observation.source_query_start
            hi = lo + len(observation.sequence)
            assert record.query_name == observation.read_id
            assert not original['source_reverse_complement']
            assert record.query_sequence[lo:hi] == observation.sequence
            assert min(record.query_qualities[lo:hi]) == original["minimum_base_quality"] >= 10
            assert not record.flag & (4 | 256 | 512 | 1024 | 2048)
            assert record.has_tag("SA")
        checked.append((entry["event"], entry["sample"]))
    assert len(checked) == 5
