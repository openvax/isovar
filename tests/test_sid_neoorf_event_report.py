"""Coordinate/provenance tests for the DNA → RNA → ORF report."""

import pytest

from examples.report_sid_neoorf_events import anatomy, clip_blocks, crossed_joins, dna_calls, literal_edit, group_key, rna_edit_description


def test_vcf_padding_distinguishes_deleted_base_and_interbase_insertion():
    deletion = literal_edit(["chr8", 100213209, "CG", "C"])
    assert (deletion["start"], deletion["end"], deletion["ref"], deletion["alt"]) == (
        100213209, 100213210, "G", "")
    insertion = literal_edit(["chr14", 55627965, "G", "GTT"])
    assert (insertion["start"], insertion["end"], insertion["alt"]) == (55627965, 55627965, "TT")
    assert literal_edit(["chr7", 100953130, "A", "dup"]) is None


@pytest.mark.parametrize("strand, expected", [("+", [22, 25]), ("-", [25, 28])])
def test_exon_number_and_cdna_offsets_follow_transcript_direction(strand, expected):
    model = dict(transcript_id="t", contig="chr1", strand=strand,
                 exons=[[100, 110], [200, 210], [300, 330]], cds_start=15, cds_end=40)
    located = anatomy(model, 302, 305)
    hit = located["exons"][0]
    assert hit["exon"] == (3 if strand == "+" else 1)
    assert hit["cdna_interval"] == expected
    assert hit["features"] == ["CDS"]
    mixed = anatomy(model, 105, 205)
    assert mixed["location"] == "exon_and_intron"
    assert mixed["overlapping_introns"] == [[110, 200]]
    boundary = anatomy(model, 110, 110)
    assert boundary["exons"][0]["at_exon_boundary"]
    assert anatomy(model, 150, 150)["location"] == "intron"


@pytest.mark.parametrize("strand, expected", [("+", [104, 108]), ("-", [112, 116])])
def test_orf_blocks_clip_on_both_genomic_strands(strand, expected):
    blocks = [dict(contig="chr1", strand=strand, query_start=10, query_end=30,
                   reference_start=100, reference_end=120)]
    clipped = clip_blocks(blocks, [14, 18])
    assert [clipped[0]["reference_start"], clipped[0]["reference_end"]] == expected
    assert clip_blocks(blocks, [30, 40]) == []


def test_dna_calls_require_dna_provenance_and_do_not_duplicate_aliases():
    vcf = dict(source_url="s", sample="T2", chrom="chr6", pos=100, vcf_id="x")
    dna = dict(source="ESVEE_PURPLE_LINX_PASS", vcf=vcf)
    rna = dict(source="CTAT_LR_RNA", vcf=dict(vcf, vcf_id="rna"))
    nominations = [dict(original=dict(original_calls=[dna, rna])),
                   dict(original=dict(original_calls=[dna]))]
    assert dna_calls(nominations) == [dna]
    assert dna_calls([dict(original=dict(original_calls=[rna]))]) == []


def test_sequence_groups_preserve_event_and_synonymous_rna_placement():
    row = dict(kind="SV_RNA_ORF", target=["alias"], geometry_id="a", amino_acids="MA",
               nucleotide_sequence="ATGGCT", ends_with_stop_codon=False)
    assert group_key(row) != group_key(dict(row, nucleotide_sequence="ATGGCC"))
    assert group_key(row) != group_key(dict(row, geometry_id="b"))
    assert group_key(row) == group_key(dict(row, source_id="another_product"))


def test_crossed_joins_use_original_path_index_not_filtered_list_position():
    join = dict(junction_index=3)
    assert crossed_joins(dict(junctions=[join], crossed_junctions=[3])) == [join]
    with pytest.raises(ValueError, match="missing"):
        crossed_joins(dict(junctions=[join], crossed_junctions=[2]))


def test_intronic_deletion_is_not_presented_as_a_mature_rna_deletion():
    event = dict(literal_edit=dict(ref="AAA", alt="", length_delta=-3),
                 anatomy=[dict(location="exon_and_intron", strand="+", exons=[])], outcomes=[])
    description = rna_edit_description(event)
    assert "mature RNA edit is unknown" in description
    assert "`AAA`" not in description
    event["anatomy"] = [dict(location="exon", strand="-", exons=[dict(at_exon_boundary=False)])]
    assert "`TTT` → `∅`" in rna_edit_description(event)


def test_published_index_counts_unresolved_events(tmp_path):
    from examples.report_sid_neoorf_events import publish

    event = dict(id="unknown", kind="SV", names=["Unknown partner"],
                 dna_origin="unknown_RNA_only_nomination", breakends=[], nominations=[],
                 original_dna_calls=[], anatomy=[], outcomes=[])
    publish(dict(events={event["id"]: event}, candidates=[], sequences={}), tmp_path)
    index = (tmp_path / "README.md").read_text()
    assert "All 1 nominations/geometry entries" in index
    assert "Unknown partner" in index
    assert "All 41" not in index


@pytest.mark.parametrize("status", ["not_assessable", "reconstruction_timeout", "no_candidate_paths"])
def test_report_preserves_unresolved_counts_and_completed_negative(tmp_path, monkeypatch, status):
    from types import SimpleNamespace
    import pyensembl
    from examples.report_sid_neoorf_events import build, publish, PRODUCTS
    from examples.sid_sv_audit.inventory import digest, write_json

    sid = next(iter(PRODUCTS))
    audit = tmp_path / "audit"
    inventory = dict(sources={sid: dict(id=sid)},
                     nominations={"n": dict(original=dict(name="Unresolved join"))},
                     geometries={"g": dict(nominations=["n"], breakends=[
                         dict(contig="chr1", position=100, retained_side="left"),
                         dict(contig="chr1", position=200, retained_side="right")])})
    write_json(audit / "inventory.json.gz", inventory)
    write_json(audit / "references.json.gz", dict(
        inventory_sha256=digest(audit / "inventory.json.gz"), models={}, unsupported_cds={},
        assignments={"g": dict(transcripts=[], excluded_transcripts=[])}))
    result = dict(request=dict(source_id=sid, geometry_id="g", orientation="forward"),
                  status=status, acquisition_status="acquired", input_limitations=[])
    if status == "no_candidate_paths":
        result.update(candidates=[], records=0, paths=0, event_paths=0, discovery={}, support_acquisition={})
    elif status == "not_assessable":
        result.update(acquisition_status="acquisition_error", reason="acquisition_error")
    else:
        result["error"] = "Per-view wall-time limit exceeded"
    result_path = audit / "results" / sid / "g.forward.json.gz"
    write_json(result_path, result)
    gtf = tmp_path / "annotation.gtf"
    gtf.write_text("# no overlapping transcripts in this fixture\n")
    genome = SimpleNamespace(gtf_path=str(gtf), transcript_fasta_paths=[],
                             gene_ids_at_locus=lambda *args: [])
    monkeypatch.setattr(pyensembl, "EnsemblRelease", lambda release: genome)
    screen = tmp_path / "screen.json.gz"
    write_json(screen, dict(scope="selected_regions", candidates=[], small_variant_outcomes=[], comparison={},
                            input_pins={str(result_path): digest(result_path), str(gtf): digest(gtf)}))
    checks = tmp_path / "checks.json"
    write_json(checks, {})
    ledger = build(screen, [audit], checks)
    outcome, = ledger["events"]["g"]["outcomes"]
    expected = 0 if status == "no_candidate_paths" else None
    assert outcome["candidates"] == outcome["records"] == outcome["paths"] == outcome["event_paths"] == expected
    assert outcome["status"] == status
    if status == "not_assessable":
        assert outcome["reason"] == "acquisition_error"
    elif status == "reconstruction_timeout":
        assert outcome["error"] == result["error"]
    publish(ledger, tmp_path / "report")
    sheet = (tmp_path / "report/events/g.md").read_text()
    assert "Unresolved join" in sheet and status in sheet
    if expected is None:
        assert "unknown / unknown / unknown" in sheet
