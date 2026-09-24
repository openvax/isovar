"""Real RNA coding hypotheses keep their unresolved annotation alternatives."""
import gzip
from hashlib import sha256
import json
from pathlib import Path

import pysam
import pytest

from isovar.fusion import fusion_from_dict, reconstruct_fusion
from tests.testing_helpers import fusion_input
from tests.data.osteosarc.expansion.references import translate

CORPUS = Path(__file__).parent / "data/fusions/coding-corpus"
ENTRIES = json.loads((CORPUS / "manifest.json").read_text())


@pytest.mark.parametrize("entry", ENTRIES, ids=lambda e: e["event"])
def test_original_coding_fusion_hypotheses(entry):
    raw = (CORPUS / entry["input"]).read_bytes()
    assert sha256(raw).hexdigest() == entry["sha256"]
    data = json.loads(gzip.decompress(raw))
    fusion, references, observations = fusion_from_dict(fusion_input(data))
    result = reconstruct_fusion(fusion, references, observations)
    assert result["status"] == entry["expected_status"] == "ambiguous"
    assert result["evidence"]["direct_fragments"] == entry["fragments"] == 2
    assert len(result["paths"][0]["translations"]) == 1  # One coding hypothesis, plus noncoding/other-CDS alternatives.
    protein, = result["paths"][0]["translations"]
    assert (protein["amino_acids"], protein["ends_with_stop_codon"]) == translate(
        fusion.sequence[protein["translation_start"]:])
    if entry["event"] == "ATP5MG--KMT2A":
        assert protein["translation_start"] == protein["cds_start"] == 8
        assert protein["amino_acids"] == "MAQFVRNLVEKTPALVNG"
        assert protein["ends_with_stop_codon"] and protein["complete_5prime"]
        assert protein["junction_in_translated_cds"] == [52, 52]
        assert protein["transcript_ids"] == ["ENST00000300688", "ENST00000524422"]
        assert len(protein["candidate_peptides"]) == 4
        assert all(r.startswith("donor_CDS_unavailable:") for r in result["reasons"])
    else:
        assert protein["translation_start"] == 2 and protein["cds_start"] is None
        assert protein["amino_acids"] == "LQTVHSIPLTINKEDDESPGLYGFLNVIVHSATGFKQSSKALQRPVASDFEPQGLSEAARWNSKENLLAGPSENDPNLF"
        assert not protein["complete_5prime"]
        assert result["reasons"] == ["junction_after_donor_CDS:ENST00000398512"]
    for peptide in protein["candidate_peptides"]:
        a, b = peptide["protein_interval"]
        assert protein["amino_acids"][a:b] == peptide["sequence"]
        assert any(3 * a < j < 3 * b for j in protein["junction_in_translated_cds"])

    header = pysam.AlignmentHeader.from_references(["chr" + str(i) for i in range(1, 23)], [300000000] * 22)
    for observation, original in zip(observations, data["original_records"]):
        read = pysam.AlignedSegment.fromstring(original["sam"], header)
        lo = observation.source_query_start
        hi = lo + len(observation.sequence)
        assert read.query_sequence[lo:hi] == fusion.sequence == observation.sequence
        assert min(read.query_qualities[lo:hi]) == original["minimum_base_quality"] >= 20
        assert not read.flag & (4 | 256 | 512 | 1024 | 2048)
        assert read.query_name == observation.read_id
