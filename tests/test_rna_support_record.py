"""Every output reports RNA support with the same record: reads, fragments, UMIs and cells."""

import pysam
import pytest

from isovar import (
    ReadCollector, export_protein_hypotheses, export_sv_rna_orfs, reconcile_allele_interpretations,
    reconstruct_fusion, run_isovar,
)
from isovar.cell_evidence import CellUmiAlleles
from isovar.cell_umi import missing_label_row
from isovar.rna_evidence import CELL_UMI_FIELDS, evidence_set, rna_support

RECORD = {"reads", "fragments", *CELL_UMI_FIELDS}


def supports(value, key=None):
    """Every support record inside an output, found by its key."""
    if isinstance(value, dict):
        if key is not None and (key.endswith("support") or key in ("ref", "alt", "other", "total", "informative")) \
                and "reads" in value:
            yield key, value
        for k, v in value.items():
            yield from supports(v, k)
    elif isinstance(value, list):
        for v in value:
            yield from supports(v, key)


def keys(value):
    if isinstance(value, dict):
        for k, v in value.items():
            yield k
            yield from keys(v)
    elif isinstance(value, list):
        for v in value:
            yield from keys(v)


def sv_result(tmp_path, labelled=False):
    from tests.test_sv_rna import Scenario, write_bam
    scenario = Scenario()
    reads = scenario.reads()
    if labelled:
        for i, r in enumerate(reads):
            r.set_tag("CB", "C%d" % (i % 3))
            r.set_tag("UB", "U%d" % i)
    return scenario.run(write_bam(tmp_path / "rna.bam", reads))


def fusion_result():
    from tests.test_fusion import example
    return reconstruct_fusion(*example(insert="GCT"))


def protein_export():
    from tests.test_protein_hypotheses import synthetic_result
    return export_protein_hypotheses([synthetic_result()], sample_id="s", source="rna.bam")


def labelled_protein_export(tmp_path, reads=None):
    """Cell/UMI-labelled export of single-cell reads; no read has the other allele."""
    from tests.test_cell_evidence import ALT, HEADER, REF, VARIANT
    from tests.test_sv_rna import write_bam
    path = write_bam(tmp_path / "cells.bam", ALT + REF if reads is None else reads, header=HEADER)
    with pysam.AlignmentFile(str(path)) as bam:
        return export_protein_hypotheses(run_isovar([VARIANT], bam), sample_id="s", source="cells.bam",
                                         cell_umi_alignment_file=bam)


def interpretations(tmp_path):
    from tests.test_allele_interpretations import INTERPRETATIONS, synthetic_bam
    with pysam.AlignmentFile(str(synthetic_bam(tmp_path))) as bam:
        return reconcile_allele_interpretations(INTERPRETATIONS, bam, sample_id="s", source="rna.bam",
                                                cell_umi_labels=True)


# Each output, and whether its cell/UMI labels were assessed.
OUTPUTS = {
    "sv_rna": (lambda tmp_path: sv_result(tmp_path, labelled=True), True),
    "sv_rna_orfs": (lambda tmp_path: export_sv_rna_orfs(sv_result(tmp_path)), True),
    "fusion": (lambda tmp_path: fusion_result(), False),
    "protein_hypotheses": (lambda tmp_path: protein_export(), False),
    "labelled_protein_hypotheses": (labelled_protein_export, True),
    "allele_interpretations": (interpretations, True),
}


@pytest.mark.parametrize("name", OUTPUTS)
def test_every_support_is_the_same_record(name, tmp_path):
    make, labelled = OUTPUTS[name]
    output = make(tmp_path)
    found = list(supports(output))
    assert found, name
    for key, support in found:
        assert RECORD <= set(support), (name, key, sorted(RECORD - set(support)))
        assert support["fragments"] <= support["reads"]
        # Labels are assessed for every support in a labelled output, including those without reads.
        assert all((support[field] is not None) == labelled for field in CELL_UMI_FIELDS), (name, key, support)
        if labelled:
            assert sum(support["label_statuses"].values()) == support["reads"]
            assert max(support["unlabeled_reads"], support["unknown_library_reads"]) <= support["reads"]
            if not support["reads"]:
                assert (support["umis"], support["cells"]) == (0, 0)
                assert not support["umis_complete"] and not support["cells_complete"]
    # Read counts are "reads" everywhere; "segments" survives only as an assembly setting.
    assert {k for k in keys(output) if "segment" in str(k)} <= {"max_extension_segments"}


def test_an_allele_without_reads_has_an_empty_evidence_set_and_zero_labels(tmp_path):
    from tests.test_cell_evidence import REF
    event, = labelled_protein_export(tmp_path, reads=REF)["events"]
    alt = event["allele_support"]["alt"]
    assert (alt["reads"], alt["umis"], alt["cells"], alt["umis_complete"]) == (0, 0, 0, False)
    assert alt["evidence_set_id"] == evidence_set(["s", "cells.bam"], [])["evidence_set_id"]
    # Only the alternate allele stores an evidence set; the others give counts only.
    assert all("evidence_set_id" not in event["allele_support"][a] for a in ("ref", "other", "total"))


def test_cell_umi_fields_are_null_only_when_not_assessed(tmp_path):
    record = rna_support([("a", "r1", 64), ("a", "r1", 128), ("a", "r2", 0)])
    assert (record["reads"], record["fragments"]) == (3, 2)
    assert all(record[field] is None for field in CELL_UMI_FIELDS)
    from tests.test_cell_evidence import ALT, HEADER, VARIANT
    from tests.test_sv_rna import write_bam
    path = write_bam(tmp_path / "cells.bam", ALT, header=HEADER)
    with pysam.AlignmentFile(str(path)) as bam:
        read_evidence = ReadCollector().read_evidence_for_variant(VARIANT, bam)
        alt = CellUmiAlleles(bam, sample_id="s", source="cells.bam").evidence(VARIANT, read_evidence)["alleles"]["alt"]
    assert RECORD <= set(alt) and all(alt[field] is not None for field in CELL_UMI_FIELDS)
    assert (alt["reads"], alt["fragments"], alt["umis"], alt["cells"]) == (4, 4, 3, 2)


def test_a_repeated_read_is_labelled_once():
    record = rna_support([("a", "r1", 0), ("a", "r1", 0)], missing_label_row)
    assert (record["reads"], record["unlabeled_reads"], record["label_statuses"]) == (
        1, 1, {"metadata_unavailable": 1})
