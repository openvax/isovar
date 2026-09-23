"""Unannotated SVs retain sequence evidence without inventing coding anchors."""

import pysam
import pytest

from isovar import export_sv_rna_orfs, reconstruct_sv_rna
from isovar.fusion import FusionBreakpoint
from tests.test_sv_rna import aligned, write_bam


@pytest.mark.parametrize("acceptor_contig", ["1", "2"])
def test_intergenic_orf_reconstruction_and_export(tmp_path, acceptor_contig):
    sequence = "CCC" + "ATG" + "GCT" * 16 + "TAA" + "CCC"
    positions = [("1", 1000 + i, "+") for i in range(30)] + [
        (acceptor_contig, 2000 + i, "+") for i in range(len(sequence) - 30)]
    records = [read for name in ("a", "b") for read in aligned(name, sequence, positions)]
    bam_path = write_bam(tmp_path / "intergenic.bam", records)
    kwargs = dict(event_id="intergenic-SV", reference_name="synthetic", references=[],
                  donor=FusionBreakpoint("1", 1030, "+"),
                  acceptor=FusionBreakpoint(acceptor_contig, 2000, "+"), regions=[],
                  sample_id="sample", source="original-RNA", event_provenance={"caller": "synthetic"},
                  assemble=False)
    with pysam.AlignmentFile(bam_path) as bam:
        result = reconstruct_sv_rna(bam, **kwargs)
    assert result["reference_models"] == []
    assert "reference_models_unavailable" in result["limitations"]
    assert result["paths"]
    assert all(p["frame_status"] == "no_gene_anchor" and not p["translations"]
               for p in result["paths"])
    candidate, = export_sv_rna_orfs(result)["candidates"]
    assert candidate["amino_acids"] == "M" + "A" * 16
    assert candidate["nucleotide_sequence"] == sequence[3:-3]
    assert candidate["rna_support"]["fragments"] == 2
    assert candidate["ends_with_stop_codon"]
    assert not candidate["initiation_observed"]
    assert not candidate["translation_observed"]
    assert not candidate["peptide_novelty_assessed"]
    assert "reconstruction_limit:reference_models_unavailable" in candidate["uncertainty_flags"]
    # Event linkage describes compatibility with the supplied adjacency;
    # without models it cannot exclude an annotated-splice explanation.
    assert result["status"] == "event_linked_candidates"
    empty_path = write_bam(tmp_path / "empty.bam", [])
    with pysam.AlignmentFile(empty_path) as bam:
        empty = reconstruct_sv_rna(bam, **kwargs)
    assert empty["status"] == "no_candidate_paths"
    assert export_sv_rna_orfs(empty)["candidates"] == []
    with pysam.AlignmentFile(bam_path) as bam, pytest.raises(ValueError, match="reference models"):
        reconstruct_sv_rna(bam, **dict(kwargs, references=[object()]))
