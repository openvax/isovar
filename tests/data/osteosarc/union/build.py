"""Rebuild the full-panel ledger from pinned packaged reads and references.

Run from the checkout with ``python -m tests.data.osteosarc.union.build``.
Use --reference-source with the original Ensembl 87 archive directory to
rebuild the reference subset first (the reference output must not exist).
No read acquisition, reselection, network access or threshold adjustment.
"""

import argparse
from collections import Counter
from hashlib import sha256
import json
import logging
from pathlib import Path
import tempfile

from osteosarc import bundle_file
import pysam
from varcode import Variant

from isovar import sid_data
from tests.data.osteosarc.expansion.references import (
    apply_variant, build_reference, load_reference, reference_genome,
)
from tests.data.osteosarc.expansion.runner import audit_mode
from tests.osteosarc_union_helpers import canonical_audit, direct_evidence, target_record


DATA = Path(__file__).parent


def build(reference_source=None):
    folder = sid_data.bundle()
    recipe = json.loads((folder / "recipe.json").read_text())
    targets = {n: t for n, t in recipe["targets"].items()
               if t["kind"] == "small_variant" or (t["kind"] == "unresolved" and t.get("label") == "current")}
    members = {v["target"]: n for n, v in recipe["members"].items() if v["source"] == sid_data.PANEL_SOURCE
               and v["target"] in targets}
    if reference_source:
        build_reference(reference_source, DATA / "reference",
                        [target_record(n, targets[n]) for n in sorted(members)])
    reference, models = load_reference(DATA / "reference")
    cases = {}
    logging.disable(logging.CRITICAL)
    with tempfile.TemporaryDirectory(prefix="isovar-union-") as work:
        genome = reference_genome(DATA / "reference", Path(work) / "reference")
        for name, target in targets.items():
            if name not in members:
                cases[name] = dict(status="unresolved_allele" if target["kind"] == "unresolved"
                                   else "no_T2_member_for_reference")
                continue
            record = target_record(name, target)
            path = bundle_file(folder, members[name], format="bam", cache=work)
            with pysam.AlignmentFile(path) as bam:
                reads = list(bam)
            direct = direct_evidence(reads, record)
            variant = Variant(record["chrom"].removeprefix("chr").replace("MT", "M"),
                              record["pos"], record["ref"], record["alt"], ensembl=genome)
            expected = {tid: apply_variant(record, models[tid]) for tid in reference["variant_transcripts"][name]}
            result = canonical_audit(audit_mode(path, variant, expected))
            if result["status"] != "ok" or any(p["validation_status"] != "ok" for p in result["proteins"]):
                raise ValueError("Cannot freeze failed pipeline/validation for %s: %s" % (name, result))
            cases[name] = dict(status="tested", member=members[name], records=len(reads),
                               direct=direct, pipeline=result)
            print(name, result.get("outcome", result["status"]), flush=True)
    manifest = dict(schema_version=1, shared_bundle=sid_data.BUNDLE,
                    shared_manifest_sha256=sid_data.BUNDLE_MANIFEST_SHA256,
                    reference_manifest_sha256=sha256((DATA / "reference/manifest.json").read_bytes()).hexdigest(),
                    source=sid_data.PANEL_SOURCE, targets=targets, cases=cases)
    (DATA / "manifest.json").write_text(json.dumps(manifest, sort_keys=True, indent=2) + "\n")
    print(Counter(c.get("pipeline", {}).get("outcome", c["status"]) for c in cases.values()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-source", type=Path)
    build(parser.parse_args().reference_source)
