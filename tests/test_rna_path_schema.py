"""Supplied fusions and reconstructed SV paths share one RNA path structure."""

import pytest

from isovar import reconstruct_fusion
from tests.test_fusion import example
from tests.test_sv_rna import Scenario, write_bam

TOP_LEVEL = {"schema", "event_id", "reference_name", "sample_id", "status", "paths", "parameters"}
PATH = {"path_id", "sequence", "junctions", "frame_status", "translations"}
JUNCTION = {"query_interval", "left", "right", "unplaced_bases", "relation", "direct_fragments"}
TRANSLATION = {"translation_start", "translation_end", "amino_acids", "ends_with_stop_codon",
               "complete_5prime", "transcript_ids", "frame_evidence", "candidate_peptides",
               "translation_observed"}


@pytest.fixture(params=["fusion", "sv-rna"])
def result(request, tmp_path):
    if request.param == "fusion":
        fusion, references, reads = example(insert="GCT")
        return reconstruct_fusion(fusion, references, reads)
    scenario = Scenario()
    return scenario.run(write_bam(tmp_path / "rna.bam", scenario.reads()))


def test_both_producers_share_the_path_junction_and_translation_fields(result):
    assert TOP_LEVEL <= set(result)
    assert isinstance(result["schema"], str)
    path = result["paths"][0]
    assert PATH <= set(path) and path["frame_status"] == "translated"
    for junction in path["junctions"]:
        assert JUNCTION <= set(junction)
    translation, = path["translations"]
    assert TRANSLATION <= set(translation) and translation["translation_observed"] is False
    start, end = translation["translation_start"], translation["translation_end"]
    assert end - start == 3 * (len(translation["amino_acids"]) + translation["ends_with_stop_codon"])
    for peptide in translation["candidate_peptides"]:
        a, b = peptide["protein_interval"]
        assert translation["amino_acids"][a:b] == peptide["sequence"]


def test_junction_intervals_are_half_open_unplaced_bases(result):
    for path in result["paths"]:
        for junction in path["junctions"]:
            start, end = junction["query_interval"]
            assert 0 < start <= end < len(path["sequence"])
            assert path["sequence"][start:end] == junction["unplaced_bases"]
