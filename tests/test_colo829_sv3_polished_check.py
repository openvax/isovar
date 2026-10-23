"""A real polished-sequence decoy must not masquerade as an SV3 join."""

import json
from pathlib import Path

from Bio import SeqIO

from reports.hcc1395_colo829.build_colo829_sv3_polished_check import _motifs, scan_sequences


DATA = Path(__file__).parent / "data/colo829_sv3"
ASSESSMENT = Path(__file__).parents[1] / "reports/hcc1395_colo829/colo829_sv3_candidate_orfs.assessment.json"


def _write_fastq(path, records):
    with path.open("w") as stream:
        for name, sequence in records:
            stream.write(f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}\n")


def test_reverse_complement_motif_is_not_full_junction_evidence(tmp_path):
    assessment = json.loads(ASSESSMENT.read_text())
    motifs = _motifs(assessment)
    decoys = [(record.id, str(record.seq)) for record in
              SeqIO.parse(DATA / "polished_reverse_motif_decoys.fasta", "fasta")]
    assert len(decoys) == 2
    fastq = tmp_path / "polished.fq"
    _write_fastq(fastq, decoys)
    count, hits = scan_sequences(fastq, motifs)
    assert count == 2
    assert hits["junction_reverse_complement"] == ["PB.1153.13", "PB.1153.23"]
    assert hits["junction_plus"] == []
    assert all(not names for key, names in hits.items() if key.startswith("window_"))

    # A synthetic positive carries the entire 19-nt join and one 27-nt
    # codon window. It must be classified separately from the two decoys.
    positive = motifs["window_plus_QAYSRTSFG"]
    _write_fastq(fastq, [*decoys, ("synthetic_true_join", positive)])
    count, hits = scan_sequences(fastq, motifs)
    assert count == 3
    assert hits["junction_plus"] == ["synthetic_true_join"]
    assert hits["window_plus_QAYSRTSFG"] == ["synthetic_true_join"]
