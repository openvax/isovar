"""Reproduce TPST1's candidate upstream ORFs from the original Sid reads.

Run from a repository checkout: python -m examples.osteosarc_sv_orfs output.json
Only reference/event metadata comes from tests/data; the RNA is exported from
osteosarc's openvax-v1 bundle (downloaded once, then offline).
"""

import argparse
import gzip
import json
from pathlib import Path
import tempfile

import pysam

from isovar import __version__
from isovar.sid_data import export
from isovar.sv_rna import reconstruct_sv_rna, sv_rna_input_from_dict


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    directory = root / "tests/data/fusions"
    entries = json.loads((directory / "long-read/manifest.json").read_text())
    entry, = [e for e in entries if (e["event"], e["source"]) == ("TPST1--CRCP", "PacBio-T1")]
    with gzip.open(directory / entry["input"], "rt") as handle:
        data = json.load(handle)
    fusion = data["fusion"]
    inputs = sv_rna_input_from_dict(dict(
        event_id=entry["event"], reference_name="GRCh38", sample_id="Sid-T1",
        donor=fusion["donor"], acceptor=fusion["acceptor"], references=data["references"],
        event_provenance=dict(fixture=entry["file"])))
    with tempfile.TemporaryDirectory(prefix="isovar-sv-orfs-") as temporary:
        sam, bam = Path(temporary) / "selected.sam", Path(temporary) / "selected.bam"
        export("fusions/long-read/" + entry["file"], sam)
        pysam.sort("-o", str(bam), str(sam))
        pysam.index(str(bam))
        with pysam.AlignmentFile(str(bam)) as alignments:
            result = reconstruct_sv_rna(alignments, source=entry["url"], **inputs)
    with args.output.open("x") as handle:
        json.dump(dict(isovar_version=__version__, result=result), handle, indent=2)
        handle.write("\n")


if __name__ == "__main__":
    main()
