"""A read set shared by several hypotheses is stored, and hashed, once."""

from isovar import rna_evidence
from isovar.allele_read import AlleleRead
from isovar.rna_evidence import EvidenceSets, evidence_set


def _read(name, placement=0):
    return AlleleRead("AC", "G", "TT", name, source_alignments=[(("group", name, 0), (0, placement, "5M", False))])


def test_repeated_read_set_is_built_once_and_matches_a_fresh_build(monkeypatch):
    built = []

    def counted(scope, identities):
        built.append(frozenset(identities))
        return evidence_set(scope, identities)

    monkeypatch.setattr(rna_evidence, "evidence_set", counted)
    a, b, c = _read("a"), _read("b"), _read("c")
    evidence = EvidenceSets(["sample", "library"])
    first = evidence.support([a, b])
    # The same segments in another order, and seen through another placement, are the same set.
    assert evidence.support([b, _read("a", placement=7), a]) == first
    other = evidence.support([a, b, c])
    assert evidence.support([]) == evidence.support([])
    assert len(built) == 3
    assert other["evidence_set_id"] != first["evidence_set_id"]

    fresh = {}
    for reads in ([a, b], [a, b, c], []):
        expected = evidence_set(["sample", "library"], [key for read in reads for key, _ in read.source_alignments])
        fresh[expected["evidence_set_id"]] = expected
    assert evidence.sets == fresh
    assert list(evidence.sets) == list(fresh)


def test_support_records_of_a_shared_set_are_independent():
    evidence = EvidenceSets(["sample", "library"])
    first = evidence.support([_read("a")])
    first["reads"] = 99
    assert evidence.support([_read("a")])["reads"] == 1


def test_reads_without_identities_store_no_set():
    evidence = EvidenceSets(["sample", "library"])
    record = evidence.support([AlleleRead("AC", "G", "TT", "legacy")])
    assert record["evidence_set_id"] is None and not evidence.sets
