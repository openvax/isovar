"""A published canine somatic call and original, region-extracted tumor RNA."""

from examples.wlos0001_tp53 import run_example


def test_wlos0001_tp53_rna_supports_mutant_protein(tmp_path):
    result = run_example(scratch_dir=tmp_path)
    event, = result["events"]
    assert event["variant"] == {
        "reference_name": "CanFam3.1", "contig": "5", "start": 32563056,
        "ref": "G", "alt": "A", "description": "chr5 g.32563056G>A",
    }
    support = event["allele_support"]
    assert [support[allele]["reads"] for allele in ("ref", "alt", "other")] == [107, 297, 1]
    assert [support[allele]["fragments"] for allele in ("ref", "alt", "other")] == [102, 289, 1]
    assert event["protein_hypotheses_complete"] is True
    supported = [
        hypothesis for hypothesis in event["protein_hypotheses"]
        if hypothesis["mutant_amino_acids"] == "C"
        and hypothesis["amino_acids"][slice(*hypothesis["mutation_interval"])] == "C"
        and hypothesis["rna_support"]["fragments"] >= 200
    ]
    assert supported
    assert any("ENSCAFT00000026465" in hypothesis["transcript_ids"] for hypothesis in supported)
