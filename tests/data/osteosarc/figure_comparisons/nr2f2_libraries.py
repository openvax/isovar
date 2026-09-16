"""Bounded independent-product NR2F2 audit; never fetch whole remote BAMs."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import gzip
from hashlib import sha256
import json
from pathlib import Path
import subprocess
from urllib.parse import quote

import pysam

from tests.data.osteosarc.expansion.acquire import assembly_from_header
from tests.data.osteosarc.expansion.inventory import fetch_snapshot
from tests.data.osteosarc.expansion.liftover import (
    chain_blocks,
    map_interval,
    genomic_sequence,
)
from tests.data.osteosarc.figure_comparisons import nr2f2


def acquire(output, source_metadata, index_cache=None):
    """Audit complete indexed regions; processing products are not replicates."""
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    repo = nr2f2.ROOT
    cache = Path(index_cache) if index_cache else None
    metadata = Path(source_metadata)
    with metadata.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    matrix = json.loads(
        gzip.decompress(
            (repo / "tests/data/osteosarc/expansion/audit/matrix.json.gz").read_bytes()
        )
    )
    inventory = {s["key"]: s for s in matrix["sources"]}
    sources = []
    for row in rows:
        if row["tissue"] == "Tumor" and "RNA" in row["assay"]:
            sources.append(dict(inventory[row["s3_path"]], sample_metadata=row))
    extra_ids = {
        "0066232879babe83",
        "53f498a544883d51",
        "1120a096937e29e1",
        "58c4d68e5f7e55b5",
        "c89442609fffd3f1",
    }
    sources += [s for s in matrix["sources"] if s["source_id"] in extra_ids]
    chain_path = out / "hg38ToHg19.over.chain.gz"
    chain_receipt = fetch_snapshot(
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz",
        chain_path,
    )
    with gzip.open(chain_path, "rt") as handle:
        blocks = list(chain_blocks(handle, {"chr15"}))
    mapping = map_interval(blocks, "chr15", 96332281, 96332371)
    assert mapping == dict(
        chrom="chr15",
        start=96875510,
        end=96875600,
        strand="+",
        chain_id=mapping["chain_id"],
    )
    original = nr2f2.load()
    assert sha256(metadata.read_bytes()).hexdigest() == original["metadata_sha256"]
    sequence, receipt = genomic_sequence(out, "hg38", "chr15", 96332281, 96332371)
    assert sequence == original["reference"]["dna"].upper()
    hg38 = json.loads((out / "hg38-chr15-96332281-96332371.json").read_text())
    references = {"GRCh37": original["reference"], "GRCh38": hg38}
    config = {
        "GRCh37": (96875527, (96875576, 96875579), (96875199, 96875850)),
        "GRCh38": (96332298, (96332347, 96332350), (96331970, 96332621)),
    }
    indices = {}
    for path in (
        sorted(cache.glob("*/inputs/alignments/*/source-index.*")) if cache else ()
    ):
        if path.suffix in (".bai", ".csi"):
            indices.setdefault(path.parent.name, []).append(path)

    def acquire_source(source):
        identity = source["source_id"]
        directory = out / identity
        directory.mkdir(exist_ok=True)
        output = directory / "evidence.json.gz"
        if output.exists():
            previous = json.loads(gzip.decompress(output.read_bytes()))
            if previous["status"] != "error":
                assert previous["key"] == source["key"]
                if previous["status"] == "ok":
                    bam_path = directory / "region.bam"
                    assert (
                        sha256(bam_path.read_bytes()).hexdigest()
                        == previous["regional_bam_sha256"]
                    )
                    with pysam.AlignmentFile(bam_path) as handle:
                        assert str(handle.header) == previous["header"]
                        assert [r.to_string() for r in handle] == previous["records"]
                return previous
        url = nr2f2.BUCKET + quote(source["key"], safe="/")
        result = dict(
            source_id=identity,
            source_url=url,
            key=source["key"],
            sample_metadata=source.get("sample_metadata"),
        )
        if not source["listed_indexes"]:
            result.update(
                status="not_assessed", reason="No index in pinned source inventory"
            )
            output.write_bytes(
                gzip.compress(json.dumps(result, sort_keys=True).encode(), mtime=0)
            )
            return result
        try:
            header_text = subprocess.check_output(
                ["samtools", "view", "--no-PG", "-H", url], timeout=90
            ).decode()
            header = pysam.AlignmentHeader.from_text(header_text)
            assembly = assembly_from_header(header.to_dict())
            assert assembly in config, assembly
            focal, deletion, (start, end) = config[assembly]
            (contig,) = [s for s in header.references if s in ("chr15", "15")]
            assert (
                header.get_reference_length(contig)
                == {"GRCh37": 102531392, "GRCh38": 101991189}[assembly]
            )
            index_url = nr2f2.BUCKET + quote(source["listed_indexes"][0], safe="/")
            index = None
            for candidate in indices.get(identity, []):
                meta = json.loads(
                    candidate.with_name(candidate.name + ".receipt.json").read_text()
                )
                if (
                    meta["url"] == index_url
                    and sha256(candidate.read_bytes()).hexdigest() == meta["sha256"]
                ):
                    index = candidate
                    index_receipt = meta
                    break
            if index is None:
                index = directory / (
                    "source-index" + Path(source["listed_indexes"][0]).suffix
                )
                index_receipt = fetch_snapshot(index_url, index)
            bam = directory / "region.bam"
            region = "%s:%d-%d" % (contig, start + 1, end)
            subprocess.run(
                [
                    "samtools",
                    "view",
                    "--no-PG",
                    "-b",
                    "-M",
                    "-X",
                    url,
                    str(index),
                    region,
                    "-o",
                    str(bam),
                ],
                check=True,
                capture_output=True,
                timeout=300,
            )
            subprocess.run(
                ["samtools", "quickcheck", str(bam)], check=True, capture_output=True
            )
            with pysam.AlignmentFile(bam) as handle:
                result.update(
                    header=str(handle.header), records=[r.to_string() for r in handle]
                )
            result.update(
                status="ok",
                assembly=assembly,
                region=region,
                index=index_receipt,
                regional_bam_sha256=sha256(bam.read_bytes()).hexdigest(),
            )
            result["audits"] = {
                name: nr2f2.recount(
                    result,
                    references[assembly],
                    min_quality=q,
                    min_mapq=m,
                    focal_position=focal,
                    deletion_interval=deletion,
                )
                for name, q, m in [
                    ("q20", 20, 20),
                    ("q30", 30, 20),
                    ("q10", 10, 20),
                    ("mapq1", 20, 1),
                ]
            }
            print(
                identity,
                source.get("sample_metadata", {}).get("display_name", source["key"]),
                len(result["records"]),
                result["audits"]["q20"]["counts"],
                flush=True,
            )
        except Exception as error:
            result.update(status="error", error=repr(error))
            if isinstance(error, subprocess.CalledProcessError):
                result["stderr"] = (
                    error.stderr.decode(errors="replace") if error.stderr else ""
                )
            print(identity, result["error"], flush=True)
        output.write_bytes(
            gzip.compress(json.dumps(result, sort_keys=True).encode(), mtime=0)
        )
        return result

    with ThreadPoolExecutor(max_workers=4) as executor:
        results = list(executor.map(acquire_source, sources))
    for source in results:
        if source["status"] == "ok":
            focal, deletion, _ = config[source["assembly"]]
            source["audits"]["sequence_only"] = nr2f2.recount(
                source,
                references[source["assembly"]],
                min_quality=None,
                focal_position=focal,
                deletion_interval=deletion,
            )
    data = dict(
        sources=results,
        references=references,
        mapping=mapping,
        chain=chain_receipt,
        reference_receipt=receipt,
        metadata_url=nr2f2.METADATA_URL,
        metadata_sha256=sha256(metadata.read_bytes()).hexdigest(),
    )
    raw = gzip.compress(json.dumps(data, sort_keys=True).encode(), mtime=0)
    (out / "nr2f2-libraries.json.gz").write_bytes(raw)
    (out / "nr2f2-libraries-manifest.json").write_text(
        json.dumps(
            dict(file="nr2f2-libraries.json.gz", sha256=sha256(raw).hexdigest()),
            indent=2,
        )
        + "\n"
    )
    errors = {s["source_id"]: s["error"] for s in results if s["status"] == "error"}
    print("COMPLETE", len(results), errors, flush=True)
    if errors:
        raise RuntimeError(errors)
    return data


def load():
    manifest = json.loads((nr2f2.CORPUS / "nr2f2-libraries-manifest.json").read_text())
    raw = (nr2f2.CORPUS / manifest["file"]).read_bytes()
    if sha256(raw).hexdigest() != manifest["sha256"]:
        raise ValueError("NR2F2 library evidence checksum mismatch")
    return json.loads(gzip.decompress(raw))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source-metadata", type=Path, required=True)
    parser.add_argument("--index-cache", type=Path)
    args = parser.parse_args()
    acquire(args.output, args.source_metadata, args.index_cache)
