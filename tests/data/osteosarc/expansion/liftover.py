"""Pin exact, ungapped UCSC chain mappings and validate original genomic alleles.

The cohort uses GRCh38 coordinates. GRCh37 nuclear coordinates are derived
from the original chain, never guessed from labels. hg19 chrM (16571 bases)
is NOT interchangeable with the 16569-base rCRS MT used by Ensembl/Tempus.
"""

import argparse
from concurrent.futures import ThreadPoolExecutor
import gzip
import json
from pathlib import Path
import sys
from urllib.parse import urlencode

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
from tests.data.osteosarc.expansion.inventory import digest, fetch_snapshot, write_json  # noqa: E402
from tests.data.osteosarc.expansion.references import reverse_complement  # noqa: E402

CHAIN_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz"


def chain_blocks(lines, contigs=None):
    """Yield forward-source, forward-destination intervals from chain blocks."""
    header = None
    for line in lines:
        fields = line.split()
        if not fields:
            continue
        if fields[0] == "chain":
            if len(fields) != 13 or fields[4] != "+" or fields[9] not in ("+", "-"):
                raise ValueError("Unsupported/malformed chain header")
            header = fields
            source, target = int(fields[5]), int(fields[10])
            continue
        if header is None or len(fields) not in (1, 3):
            raise ValueError("Chain block without valid header")
        size = int(fields[0])
        if size <= 0:
            raise ValueError("Nonpositive chain block")
        reverse = header[9] == "-"
        target_start = int(header[8]) - target - size if reverse else target
        if contigs is None or header[2] in contigs:
            yield dict(source_contig=header[2], source_start=source, source_end=source + size,
                       target_contig=header[7], target_start=target_start, target_end=target_start + size,
                       strand=header[9], chain_id=header[12])
        source += size
        target += size
        if len(fields) == 3:
            source += int(fields[1])
            target += int(fields[2])
        else:
            if source != int(header[6]) or target != int(header[11]):
                raise ValueError("Chain endpoint does not match blocks")
            header = None
    if header is not None:
        raise ValueError("Truncated chain")


def map_interval(blocks, chrom, start, end):
    """Require one whole-allele ungapped block; reject gaps and ambiguity."""
    if start < 0 or end <= start:
        raise ValueError("Invalid source interval")
    matches = []
    for block in blocks:
        if block["source_contig"] != chrom or not block["source_start"] <= start < end <= block["source_end"]:
            continue
        if block["strand"] == "+":
            a = block["target_start"] + start - block["source_start"]
            b = a + end - start
        else:
            b = block["target_end"] - (start - block["source_start"])
            a = b - (end - start)
        matches.append(dict(chrom=block["target_contig"], start=a, end=b,
                            strand=block["strand"], chain_id=block["chain_id"]))
    if len(matches) != 1:
        raise ValueError(f"Expected one ungapped interval mapping, found {len(matches)}")
    return matches[0]


def lift_variant(record, blocks):
    mapping = map_interval(blocks, record["chrom"], record["pos"] - 1, record["pos"] - 1 + len(record["ref"]))
    ref, alt = record["ref"], record["alt"]
    if mapping["strand"] == "-":
        # Reverse-strand indels need a new genomic left anchor. Do not emit
        # invalid VCF or conceal this as a successfully mapped substitution.
        if len(ref) != len(alt):
            raise ValueError("Reverse-chain indel requires independently validated reanchoring")
        ref, alt = reverse_complement(ref), reverse_complement(alt)
    return dict(record, assembly="GRCh37", chrom=mapping["chrom"], pos=mapping["start"] + 1,
                ref=ref, alt=alt, original_identity={k: record[k] for k in ("assembly", "chrom", "pos", "ref", "alt")},
                coordinate_mapping=mapping)


def genomic_sequence(directory, genome, chrom, start, end):
    """Fetch bounded original UCSC sequence with identity and length checks."""
    path = Path(directory) / f"{genome}-{chrom}-{start}-{end}.json"
    query = urlencode(dict(genome=genome, chrom=chrom, start=start, end=end))
    receipt = fetch_snapshot("https://api.genome.ucsc.edu/getData/sequence?" + query, path)
    response = json.loads(path.read_text())
    if any(response.get(k) != value for k, value in dict(genome=genome, chrom=chrom, start=start, end=end).items()):
        raise ValueError("Genomic response coordinate identity mismatch")
    sequence = response["dna"].upper()
    if len(sequence) != end - start:
        raise ValueError("Genomic response sequence length mismatch")
    return sequence, receipt


def verify_allele(record, genome, directory):
    start = max(0, record["pos"] - 1 - 30)
    end = record["pos"] - 1 + len(record["ref"]) + 30
    sequence, receipt = genomic_sequence(directory, genome, record["chrom"], start, end)
    offset = record["pos"] - 1 - start
    observed = sequence[offset:offset + len(record["ref"])]
    if observed != record["ref"]:
        raise ValueError(f"Genomic reference mismatch for {record['variant_id']}: {observed} != {record['ref']}")
    return dict(record, genomic_validation=dict(genome=genome, start=start, end=end, sequence=sequence,
                                                allele_offset=offset, receipt=receipt))


def prepare(destination):
    destination = Path(destination)
    inventory = json.loads((destination / "inventory-live.json").read_text())
    directory = destination / "genomic-references"
    directory.mkdir(exist_ok=True)
    chain = directory / "hg38ToHg19.over.chain.gz"
    receipt = fetch_snapshot(CHAIN_URL, chain)
    with gzip.open(chain, "rt") as handle:
        blocks = list(chain_blocks(handle, {v["chrom"] for v in inventory["variants"]}))
    with ThreadPoolExecutor(max_workers=4) as executor:
        original = list(executor.map(lambda r: verify_allele(r, "hg38", directory), inventory["variants"]))
        lifted = [lift_variant(r, blocks) for r in inventory["variants"]]
        hg19 = list(executor.map(lambda r: verify_allele(r, "hg19", directory), lifted))
    # Explicitly retain both mitochondrial coordinate systems. For rCRS,
    # unchanged positions are justified by the pinned hg38 MT sequence and
    # must also match the independently pinned Ensembl GRCh37 transcript.
    grch37 = [dict(o, assembly="GRCh37", original_identity={k: o[k] for k in ("assembly", "chrom", "pos", "ref", "alt")},
                   coordinate_mapping=dict(method="rCRS MT, not UCSC hg19 chrM", reference_length=16569))
              if o["chrom"] == "chrM" else n for o, n in zip(original, hg19)]
    for name, variants in (("GRCh38", original), ("GRCh37", grch37), ("GRCh37-hg19MT", hg19)):
        value = dict(inventory, variants=variants,
                     coordinate_provenance=dict(chain=receipt, source_inventory_sha256=digest(destination / "inventory-live.json")))
        path = destination / f"inventory-{name}-validated.json"
        if path.exists():
            if json.loads(path.read_text()) != value:
                raise ValueError("Validated inventory drift")
        else:
            write_json(path, value)
        print("Validated genomic alleles", name, len(variants), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", type=Path)
    prepare(parser.parse_args().destination)
