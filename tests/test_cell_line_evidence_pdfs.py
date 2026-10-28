"""Check the published cell-line report counts and reference comparisons."""

from pathlib import Path

import pytest


@pytest.mark.parametrize(
    ("sample", "sv_counts", "example_calls"),
    [
        ("HCC1395", (1788, 354, 172, 6, 3, 0), ("COSM32759", "COSM32783")),
        ("COLO829", (188, 97, 64, 6, 0, 0),
         ("chr9:21971154:CCG>C", "chr2:105361307:G>GCCCC",
          "chr1:156265990:G>A")),
    ],
)
def test_cell_line_sv_funnels_and_protein_comparisons(tmp_path, sample, sv_counts,
                                                      example_calls):
    pypdf = pytest.importorskip("pypdf")
    from reports.hcc1395_colo829.build_evidence_pdfs import build_one, collect

    data = collect(sample)
    assert tuple(data["counts"][name] for name in (
        "sv_total", "sv_exonic", "sv_exonic_covered", "sv_exonic_candidate",
        "sv_exonic_junction", "sv_exonic_product")) == sv_counts
    assert set(example_calls) <= {row['call_id'] for row in data['example_references']}
    for row in data["example_references"]:
        mutant = data["summary"]["top"][row["call_id"]]
        start = int(mutant["mutation_start_idx"])
        if row['placement_method'] == 'unique_shared_protein_prefix':
            assert row["reference_context"][:start] == mutant["protein_sequence"][:start]
        elif row['placement_method'] == 'verified_in_frame_deletion':
            left, right = int(row['deleted_reference_start']), int(row['deleted_reference_end'])
            assert row['reference_context'][:left] + row['reference_context'][right:] == mutant['protein_sequence']
        else:
            assert row['placement_method'] == 'genomic_codon_anchor_RNA_flanks_differ'
        if row['call_id'] in example_calls:
            assert row["reference_context"] != mutant["protein_sequence"]

    build_one(sample, tmp_path)
    pdf = pypdf.PdfReader(Path(tmp_path) / f"{sample}_somatic_RNA_evidence.pdf")
    cover = pdf.pages[0].extract_text()
    assert "Small variants: DNA to RNA and local protein" in cover
    assert "Structural variants: complementary routes" in cover
    assert "Soft clips / unmapped" in cover
    assert "not an eligibility gate" in cover
    full_text = "\n".join(page.extract_text() for page in pdf.pages)
    assert "RNA link candidates per exon-involved SV" in full_text
    assert "Original-read replay: junction versus RNA path" in full_text
    assert "one read need not span the assembled context" in full_text
    if sample == "COLO829":
        assert "1,213" in full_text
    assert 'Variant landscape' in full_text
    assert 'Highest-support SNVs' in full_text
    assert 'Top RNA-supported frameshift nominations' in full_text
    assert 'ncORFs and other antigen-source categories' in full_text
    assert 'IGV (GRCh38)' in full_text
    assert 'Focal change' in full_text
    if sample == 'HCC1395':
        assert 'Unclassified' in full_text or 'unclassified provenance' in full_text
        assert '531-aa' in full_text and '149 exact Illumina junction' in full_text
    else:
        assert '5,248 fragments contribute' in full_text and 'CNIH3' in full_text
        assert '1,027 fragments exactly span' in full_text
        assert 'DDX17: intronic ALT RNA' in full_text
        assert '27,309 pairs' in full_text
        assert 'chr22:38490360-38490362' in full_text
        assert 'vertical marker shows the missing reference residue' in full_text


def test_in_frame_deletion_keeps_unchanged_suffix_aligned():
    from reports.hcc1395_colo829.report_landscape import NeutralReferenceComparison
    from reports.hcc1395_colo829.build_evidence_pdfs import collect

    data = collect('COLO829')
    ref = next(r for r in data['example_references'] if r['call_id'] == 'chr14:92811641:GAGA>G')
    product = data['summary']['top'][ref['call_id']]
    assert ref['placement_method'] == 'verified_in_frame_deletion'
    deletion = int(ref['deleted_reference_start']), int(ref['deleted_reference_end'])
    comparison = NeutralReferenceComparison(ref['reference_context'], product['protein_sequence'],
                                           8, 9, deleted_reference_interval=deletion)
    assert comparison.aligned_mutant == 'LAERKYSDEK RVDELQQQVKLYKLNLESSKQEL'
    assert comparison.reference == 'LAERKYSDEKKRVDELQQQVKLYKLNLESSKQEL'
    assert int(ref['protein_start_1based']) + deletion[0] == 405


def test_deletion_gap_does_not_relabel_unexplained_changes():
    from reports.hcc1395_colo829.report_landscape import NeutralReferenceComparison

    comparison = NeutralReferenceComparison('MABCDEK', 'MADEK', 2, 2,
                                           deleted_reference_interval=(2, 4))
    assert comparison.aligned_mutant == 'MA  DEK'
    with pytest.raises(ValueError, match='preserve both observed RNA flanks'):
        NeutralReferenceComparison('MABCDEK', 'MAXEK', 2, 2,
                                   deleted_reference_interval=(2, 4))


def test_ddx17_original_rna_audit_preserves_ambiguous_junction(tmp_path):
    import json
    import shutil

    from reports.hcc1395_colo829.build_ddx17_audit import HERE, build

    directory = tmp_path / 'ddx17'
    shutil.copytree(HERE, directory)
    result = build(directory)
    assert result == json.loads((HERE / 'results.json').read_text())
    assert result['coding_reference_contexts'] == 0
    assert result['original_allele_fragments']['alt'] == 221
    assert result['counts']['alt_mapped_fragments'] == 219
    assert result['counts']['alt_locally_unique_fragments'] == 214
    assert result['fragments_in_multiple_allele_groups'] == 11
    assert result['intron_uncovered_bases'] == 0
    assert result['intron_min_locally_unique_fragment_depth'] == 115
    junctions = {(r['intron_start_0based'], r['intron_end_0based']): r
                 for r in result['junctions']}
    assert len(junctions[38488115, 38492055]['locally_unique_fragments']) == 449
    alternative = junctions[38490436, 38492055]
    assert len(alternative['locally_unique_fragments']) == 10
    assert alternative['genomic_motif'] == 'CT-AC'
    assert not alternative.get('alt_fragments')
    ambiguous = junctions[38486939, 38489622]
    assert ambiguous['alt_fragments'] == ['SRR8615617.75487984']
    assert not ambiguous.get('locally_unique_fragments')
    assert [9, 92] in ambiguous['anchor_pairs']
    assert all(h['mutation_downstream_of_stop'] for h in result['reference_frame_hypotheses'])
    assert {tuple(h['first_stop_1based']) for h in result['reference_frame_hypotheses']} == {
        (38492042, 38492044), (38490360, 38490362), (38491075, 38491077)}
    (directory / 'models.json').write_text('{}')
    with pytest.raises(ValueError, match='Input digest mismatch: models.json'):
        build(directory)
