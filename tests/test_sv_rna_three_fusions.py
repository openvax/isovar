"""Original RNA regressions with reproducible selection and annotation builders."""
import json

import pysam
import pytest

from isovar.sid_data import sam_digest
from tests.data.fusions.audit_three_fusions import DATA, audit_candidates, pinned_json, run_case
from tests.data.fusions.build_three_fusions import select
from tests.data.osteosarc.expansion.references import translate

MANIFEST = json.loads((DATA / "manifest.json").read_text())
PARD = "PARD3B--CDKN2B-AS1-CDKN2B"
GABBR = "GABBR1--SLC29A1"
OTUD = "OTUD7A--FMN1"
GABBR_PEPTIDE = "MKRLVSSSRAWWRMPVIPAPTEAEAGESLESGRRRLQ"
OTUD_PEPTIDE = "MAAKAGTSLEARSLRPAWPTC"
PARD_PEPTIDE = "MNISNIHISTQKKKKKKSRF"


@pytest.fixture(scope="module")
def case(tmp_path_factory):
    cache = {}

    def get(event, source, orientation):
        key = event, source, orientation
        if key not in cache:
            result, records = run_case(DATA, *key, tmp_path_factory.mktemp("three-fusions"))
            cache[key] = result, audit_candidates(result, records)
        return cache[key]
    return get


@pytest.mark.parametrize("entry", MANIFEST["fixtures"], ids=lambda e: e["file"])
def test_original_records_and_selection_survive_fixture_roundtrip(entry, tmp_path):
    inventory = pinned_json(DATA, MANIFEST["source_inventory"])
    assert MANIFEST["sources"][entry["source"]]["source_id"] in inventory
    data = pinned_json(DATA, entry)
    header = pysam.AlignmentHeader.from_dict(data["header"])
    unsorted, bam = tmp_path / "unsorted.bam", tmp_path / "reads.bam"
    with pysam.AlignmentFile(str(unsorted), "wb", header=header) as out:
        for digest, line in data["records"].items():
            assert sam_digest(line) == digest
            read = pysam.AlignedSegment.fromstring(line, header)
            assert read.to_string() == line  # No sequence/QUAL/tag rewriting.
            out.write(read)
    pysam.sort("--no-PG", "-o", str(bam), str(unsorted))
    pysam.index(str(bam))
    with pysam.AlignmentFile(str(bam)) as handle:
        selected, segments = select(handle, MANIFEST["events"][entry["event"]])
    assert set(selected) == set(data["records"].values())
    assert len(selected) == entry["records"] and segments == entry["segments"]


@pytest.mark.parametrize("event,source,orientation,peptide,native,labels,tail", [
    (GABBR, "ONT-T1", "forward", GABBR_PEPTIDE, 1, 1, 1),
    (GABBR, "ONT-T2", "forward", GABBR_PEPTIDE, 3, 3, 3),
    (OTUD, "ONT-T2", "reverse", OTUD_PEPTIDE, 9, 7, 8),
    (PARD, "ONT-T1", "reverse", PARD_PEPTIDE, 4, 4, 4),
    (PARD, "ONT-T3", "reverse", PARD_PEPTIDE, 2, 2, 2),
])
def test_complete_orfs_have_original_native_witnesses_but_no_proven_initiation(
        case, event, source, orientation, peptide, native, labels, tail):
    result, candidates = case(event, source, orientation)
    candidate, = [c for c in candidates if c["amino_acids"] == peptide]
    assert translate(candidate["nucleotides"], table=1) == (peptide, True)
    assert candidate["native_fragments"] == candidate["native_q10_fragments"] == native
    assert candidate["native_cb_ub_labels"] == labels
    assert candidate["native_terminal_a10_fragments"] == tail
    assert not candidate["annotated_start"]
    assert not candidate["initiation_observed"] and not candidate["translation_observed"]
    assert not any("breakpoint_junction" in t["departure_relations"]
                   for p in result["paths"] for t in p["translations"])


def test_gabbr_opposite_strand_junction_read_is_not_a_second_native_orf_witness(case):
    _, candidates = case(GABBR, "ONT-T1", "forward")
    candidate, = [c for c in candidates if c["amino_acids"] == GABBR_PEPTIDE]
    assert candidate["orientation_counts"] == {"native": 1, "opposite": 1}


def test_otud_catalogue_orientation_only_produces_reverse_complement_candidates(case):
    _, candidates = case(OTUD, "ONT-T2", "forward")
    assert candidates and any(len(c["amino_acids"]) >= 60 for c in candidates)
    assert all(c["native_fragments"] == 0 for c in candidates)
    assert all(c["orientation_counts"] == {"opposite": len(c["witnesses"])} for c in candidates)


def test_pard_repeat_alternatives_are_retained_instead_of_replaced_by_one_consensus(case):
    _, first = case(PARD, "ONT-T1", "reverse")
    _, second = case(PARD, "ONT-T2", "reverse")
    assert {c["amino_acids"]: c["native_fragments"] for c in first if c["native_fragments"]} == {
        "MNISNIHISTQKKKKKE": 4,
        "MNISNIHISTQKKKKKKKE": 3,
        PARD_PEPTIDE: 4,
        "MNISNIHISTQKKKKKRVDFKCSHHERNKYVR": 4,
        "MNISNIHISTQKKKKKSRF": 7,
    }
    assert {c["amino_acids"]: c["native_fragments"] for c in second if c["native_fragments"]} == {
        "MNISNIHISTQKKKKKE": 3,
        "MNISNIHISTQKKKKKRVDFKCSHHERNKYVR": 3,
        "MNISNIHISTQKKKKRVDFKCSHHERNKYVR": 3,
    }


def test_pard_pacbio_corroboration_preserves_missing_qualities(case):
    _, candidates = case(PARD, "PacBio-T1", "reverse")
    candidate, = [c for c in candidates if c["amino_acids"] == PARD_PEPTIDE]
    assert candidate["native_fragments"] == candidate["native_missing_quality_fragments"] == 2
    assert candidate["native_q10_fragments"] == 0
    assert candidate["native_cb_ub_labels"] == 0  # Original XM is not a UB label.
