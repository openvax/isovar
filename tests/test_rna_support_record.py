"""Every output reports RNA support with the same record: reads, fragments, UMIs and cells."""

import pysam
import pytest

from isovar import (
    ReadCollector, export_protein_hypotheses, export_sv_rna_orfs, reconcile_allele_interpretations,
    reconstruct_fusion,
)
from isovar.cell_evidence import CellUmiAlleles
from isovar.rna_evidence import CELL_UMI_FIELDS, rna_support

RECORD = {"reads", "fragments", *CELL_UMI_FIELDS}


def supports(value, key=None):
    """Every support record inside an output, found by its key."""
    if isinstance(value, dict):
        if key is not None and (key.endswith("support") or key in ("ref", "alt", "other", "total")) \
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


def sv_result(tmp_path):
    from tests.test_sv_rna import Scenario, write_bam
    scenario = Scenario()
    return scenario.run(write_bam(tmp_path / "rna.bam", scenario.reads()))


def fusion_result():
    from tests.test_fusion import example
    return reconstruct_fusion(*example(insert="GCT"))


def protein_export():
    from tests.test_protein_hypotheses import synthetic_result
    return export_protein_hypotheses([synthetic_result()], sample_id="s", source="rna.bam")


def interpretations(tmp_path):
    from tests.test_allele_interpretations import INTERPRETATIONS, synthetic_bam
    with pysam.AlignmentFile(str(synthetic_bam(tmp_path))) as bam:
        return reconcile_allele_interpretations(INTERPRETATIONS, bam, sample_id="s", source="rna.bam",
                                                cell_umi_labels=True)


OUTPUTS = {
    "sv_rna": sv_result,
    "sv_rna_orfs": lambda tmp_path: export_sv_rna_orfs(sv_result(tmp_path)),
    "fusion": lambda tmp_path: fusion_result(),
    "protein_hypotheses": lambda tmp_path: protein_export(),
    "allele_interpretations": interpretations,
}


@pytest.mark.parametrize("name", OUTPUTS)
def test_every_support_is_the_same_record(name, tmp_path):
    output = OUTPUTS[name](tmp_path)
    found = list(supports(output))
    assert found, name
    for key, support in found:
        assert RECORD <= set(support), (name, key, sorted(RECORD - set(support)))
        assert support["fragments"] <= support["reads"]
    # Read counts are "reads" everywhere; "segments" survives only as an assembly setting.
    assert {k for k in keys(output) if "segment" in str(k)} <= {"max_extension_segments"}


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
