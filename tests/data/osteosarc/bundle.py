"""Compile the reviewed test selections into an osteosarc acquisition recipe.

This inventories test inputs, not acquisition directories. The package contains
only their original reads, shared across fixtures. Report bodies, BAM indices,
reference archives and background locus downloads are not packaged.

python -m tests.data.osteosarc.bundle --snapshot NAME --cache CACHE --output RECIPE
"""

import argparse
from importlib.metadata import version
import json
from pathlib import Path
import re

import pysam

from isovar.sid_data import asset_identity, open_dataset, read_json, sam_digest, write_json
from tests.data.osteosarc.expansion.acquire import assembly_from_header


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"


def flatten_records(value):
    """Only explicit original SAM fields; never reconstructed cDNA/windows."""
    if isinstance(value, str):
        if len(value.split("\t")) < 11:
            raise ValueError("Not an original SAM record")
        yield value
    elif isinstance(value, list):
        for item in value:
            yield from flatten_records(item)
    elif isinstance(value, dict):
        for key in ("sam", "partner_sam"):
            if value.get(key):
                yield from flatten_records(value[key])
    else:
        raise ValueError("Unexpected original-record representation")


def alignment(path, url, tests):
    with pysam.AlignmentFile(DATA / path) as handle:
        return dict(name=path, url=url, records=[r.to_string() for r in handle],
                    tests=tests, header=handle.header.to_dict(), record_references=False)


def embedded(path, tests):
    value = read_json(DATA / path)
    sources = {s["id"]: s["url"] for s in value.get("sources", [])
               if isinstance(s, dict) and "id" in s and "url" in s}
    default_url = value.get("fusion", {}).get("provenance", {}).get("source")

    def walk(node, pointer="", url=default_url, header=None):
        if isinstance(node, dict):
            source = node.get("source")
            if isinstance(source, dict):
                url = source.get("url", url)
            elif source in sources:
                url = sources[source]
            url = node.get("source_url", node.get("url", node.get("receipt", {}).get("url", url)))
            header = node.get("header", header)
            fields = ("records",) if path.endswith(("nr2f2-evidence.json.gz", "nr2f2-libraries.json.gz")) else (
                "original_records", "original_sam")
            for field in fields:
                if field not in node or not node[field]:
                    continue
                # Only the reviewed inputs below use 'records' as SAM strings.
                lines = list(flatten_records(node[field]))
                if not lines:
                    continue
                if not url or "sid-sijbrandij-osteosarc-dataset" not in url:
                    raise ValueError("Missing Sid source at " + path + pointer)
                yield dict(name=path + "#" + pointer + "/" + field, url=url,
                           records=lines, tests=tests, header=header,
                           record_references=isinstance(node[field][0], (dict, list)))
            for key, item in node.items():
                if key not in ("original_records", "original_sam", "records", "header"):
                    yield from walk(item, pointer + "/" + key, url, header)
        elif isinstance(node, list):
            for index, item in enumerate(node):
                yield from walk(item, pointer + "/" + str(index), url, header)

    yield from walk(value)


def fixture_inputs():
    """Explicit, reviewable scope, with each selection's consuming tests."""
    legacy = read_json(DATA / "osteosarc/manifest.json")
    for row in legacy["datasets"].values():
        yield alignment("osteosarc/" + row["file"], row["source"]["url"],
                        ["test_osteosarc_audit.py", "test_osteosarc_proteins.py"])
    for directory, tests in (
        ("osteosarc/expansion/corpus", ["test_osteosarc_expansion_corpus.py"]),
        ("osteosarc/expansion/stress-corpus", ["test_osteosarc_expansion_stress.py"]),
        ("osteosarc/figure_comparisons/corpus", ["test_osteosarc_figures.py"]),
    ):
        for row in read_json(DATA / directory / "manifest.json")["cases"]:
            for field in ("bam", "primary_bam"):
                if field in row:
                    yield alignment(directory + "/" + row[field], row["source_url"], tests)
    pacbio = read_json(DATA / "osteosarc/figure_comparisons/corpus/pacbio-indels-manifest.json")
    yield alignment("osteosarc/figure_comparisons/corpus/" + pacbio["file"],
                    pacbio["source"]["url"], ["test_insertion_boundaries.py"])
    for row in read_json(DATA / "fusions/long-read/manifest.json"):
        yield alignment("fusions/long-read/" + row["file"], row["url"], ["test_sv_rna.py"])
    for kind in ("ont", "short"):
        urls = re.findall(r"https://\S+\.bam\b", (DATA / "chimeric/README.md").read_text())
        yield alignment("chimeric/osteosarc-" + kind + ".sam", urls[0 if kind == "ont" else 1],
                        ["test_chimeric_phasing.py"])
    urls = re.findall(r"https://[^)\s]+\.bam\b", (DATA / "read_ends/README.md").read_text())
    with pysam.AlignmentFile(DATA / "read_ends/osteosarc.sam") as handle:
        for read in handle:
            source_index = 2 if read.query_name.startswith("molecule/") else (1 if read.is_paired else 0)
            yield dict(name="read_ends/osteosarc.sam#" + read.query_name, url=urls[source_index],
                       records=[read.to_string()], header=handle.header.to_dict(),
                       tests=["test_read_end_inference.py"], record_references=False)
    for path in sorted((DATA / "fusions/corpus").glob("*.input.json.gz")):
        yield from embedded(str(path.relative_to(DATA)), ["test_fusion.py", "test_sv_rna.py"])
    yield from embedded("fusions/coding-corpus/ATP5MG--KMT2A.input.json.gz",
                        ["test_coding_fusion_examples.py", "test_sv_rna.py"])
    for name, tests in (
        ("context-evidence", ["test_aligned_reading_frame.py", "test_osteosarc_figures.py"]),
        ("rna-footprints", ["test_osteosarc_footprints.py"]),
        ("extended-footprints", ["test_extended_footprints.py"]),
        ("dlg5", ["test_extended_footprints.py"]),
        ("insertion-boundaries", ["test_insertion_boundaries.py"]),
        ("nr2f2-evidence", ["test_nr2f2_evidence.py"]),
        ("nr2f2-libraries", ["test_nr2f2_libraries.py"]),
    ):
        yield from embedded("osteosarc/figure_comparisons/corpus/" + name + ".json.gz", tests)


def compile_recipe(dataset):
    sources, fixtures, positions = {}, {}, {}
    for fixture in fixture_inputs():
        for test in fixture["tests"]:
            if not (ROOT / test).is_file():
                raise ValueError("Unknown consuming test: " + test)
        asset = dataset.file(fixture["url"])
        identity = asset.id[:16]
        source = sources.setdefault(identity, dict(asset=asset_identity(asset), assembly=None,
                                                  reference_lengths={}))
        header = fixture["header"]
        if isinstance(header, str):
            header = pysam.AlignmentHeader.from_text(header).to_dict()
        if header:
            assembly = assembly_from_header(header)
            if assembly:
                if source["assembly"] not in (None, assembly):
                    raise ValueError("Contradictory source assemblies")
                source["assembly"] = assembly
            for row in header.get("SQ", []):
                previous = source["reference_lengths"].setdefault(row["SN"], row["LN"])
                if previous != row["LN"]:
                    raise ValueError("Contradictory reference length")
        coordinates = positions.setdefault(identity, set())
        for line in fixture["records"]:
            fields = line.split("\t")
            if int(fields[1]) & 4 or fields[2] == "*" or int(fields[3]) < 1:
                raise ValueError("Unmapped input needs an explicit mate-recovery recipe: " + fixture["name"])
            coordinates.add((fields[2], int(fields[3]) - 1))
        if fixture["name"] in fixtures:
            raise ValueError("Duplicate fixture name")
        fixtures[fixture["name"]] = dict(source=identity, records=list(map(sam_digest, fixture["records"])),
                                         tests=fixture["tests"], record_references=fixture["record_references"])
    for identity, source in sources.items():
        if source["assembly"] is None:
            info = dataset.inspect_alignment(source["asset"]["key"])
            source["assembly"] = info.assembly
            source["reference_lengths"] = {s["SN"]: s["LN"] for s in info.header["SQ"]}
        if source["assembly"] not in ("GRCh37", "GRCh38"):
            raise ValueError("Unresolved source reference: " + identity)
        # One mapped base per required alignment suffices to retrieve it whole.
        # Coalesce adjacent query points, never bridge introns or broad windows.
        intervals = []
        for contig, start in sorted(positions[identity]):
            if intervals and intervals[-1][0] == contig and intervals[-1][2] == start:
                intervals[-1][2] = start + 1
            else:
                intervals.append([contig, start, start + 1])
        source["regions"] = [dict(contig=c, start=a, end=b, assembly=source["assembly"],
                                  **({"reference_length": source["reference_lengths"][c]}
                                     if c in ("M", "MT", "chrM", "chrMT") else {}))
                             for c, a, b in intervals]
        del source["reference_lengths"]
    return dict(schema_version=1, snapshot_id=dataset.id, osteosarc_version=version("osteosarc"),
                license="CC0-1.0", selection="Exact original records exercised by the listed tests; not VAF sampling",
                sources=sources, fixtures=fixtures)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snapshot", required=True)
    parser.add_argument("--cache")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    recipe = compile_recipe(open_dataset(args.snapshot, args.cache))
    write_json(args.output, recipe)
    print(json.dumps(dict(fixtures=len(recipe["fixtures"]), sources=len(recipe["sources"]),
                          references=sum(len(f["records"]) for f in recipe["fixtures"].values())), sort_keys=True))


if __name__ == "__main__":
    main()
