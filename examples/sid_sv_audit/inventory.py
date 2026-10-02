"""Freeze nominations, retained-flank geometries, RNA products and references.

The manifest is a denominator, not a positive-call filter. Original nominations
and every alignment product remain inspectable, including exclusions.
"""

from collections import Counter, defaultdict
from dataclasses import asdict
import csv
import gzip
from hashlib import sha256
import io
import json
from pathlib import Path
import re


RNA_ASSAYS = {"rna-seq", "scrna-seq", "cite-seq"}
FUSION_TABLES = ("fusion_sequences", "scan_summary")
FUSION_URL = "https://osteosarc.com/fusions/tables/"


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def identity(value):
    return sha256(canonical(value)).hexdigest()


def digest(path):
    with Path(path).open("rb") as handle:
        result = sha256()
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            result.update(chunk)
    return result.hexdigest()


def write_json(path, value):
    """Deterministic JSON; existing content may only be reused unchanged."""
    path = Path(path)
    data = canonical(value) + b"\n"
    if path.suffix == ".gz":
        data = gzip.compress(data, mtime=0)
    if path.exists():
        if path.read_bytes() != data:
            raise ValueError("Refusing to replace different audit content: %s" % path)
    else:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("xb") as handle:
            handle.write(data)


def read_json(path):
    path = Path(path)
    data = path.read_bytes()
    return json.loads(gzip.decompress(data) if path.suffix == ".gz" else data)


def geometry(breakends):
    """Canonical undirected adjacency; retained sides determine both RNA views."""
    if len(breakends) != 2:
        raise ValueError("An adjacency requires exactly two breakends")
    ends = []
    for end in breakends:
        contig, position, side = (end[k] for k in ("contig", "position", "retained_side"))
        if not contig or type(position) is not int or position < 0 or side not in ("left", "right"):
            raise ValueError("Invalid retained-flank geometry")
        ends.append(dict(contig=contig, position=position, retained_side=side))
    return sorted(ends, key=lambda end: (end["contig"], end["position"], end["retained_side"]))


def oriented_breakpoints(ends, reverse=False):
    """The donor is traversed towards its boundary; the acceptor away from it."""
    donor, acceptor = reversed(ends) if reverse else ends
    return (dict(contig=donor["contig"], position=donor["position"],
                 strand="+" if donor["retained_side"] == "left" else "-"),
            dict(contig=acceptor["contig"], position=acceptor["position"],
                 strand="+" if acceptor["retained_side"] == "right" else "-"))


def fusion_breakends(row):
    """Survey coordinates name the last retained base, one-based inclusive.

    Convert either side independently: retained-left ends after that base;
    retained-right starts before it. Do not infer strands from gene names.
    """
    ends = []
    for suffix in ("D", "A"):
        side, position = row["retain" + suffix], int(row["pos" + suffix])
        if position < 1 or side not in ("left", "right"):
            raise ValueError("Invalid fusion-survey breakend")
        ends.append(dict(contig=row["chr" + suffix],
                         position=position - (side == "right"), retained_side=side))
    return geometry(ends)


def nominations(catalogue, fusion_rows, summary_rows):
    """Keep all shared targets and every survey call; group only query geometry.

    Grouping saves reconstruction work. It does NOT declare identical alleles,
    independent mutations, or add support between callers and inserted alleles.
    """
    targets, groups = {}, {}

    def add(name, entry, ends):
        if name in targets:
            raise ValueError("Duplicate nomination identity: " + name)
        entry = dict(entry, geometry_id=None)
        if ends is not None:
            group_id = "adj-" + identity(ends)[:24]
            group = groups.setdefault(group_id, dict(breakends=ends, nominations=[]))
            if group["breakends"] != ends:
                raise ValueError("Geometry identity collision")
            group["nominations"].append(name)
            entry["geometry_id"] = group_id
        targets[name] = entry

    for name, target in sorted(catalogue["targets"].items()):
        if target["assembly"] != "GRCh38":
            raise ValueError("Unsupported shared-catalogue assembly")
        ends = None
        if target["kind"] == "sv":
            if target["coordinates"] != "zero-based-interbase":
                raise ValueError("Unrecognized shared SV coordinates")
            ends = geometry(target["breakends"])
        elif target["kind"] != "unresolved":
            raise ValueError("Unrecognized shared nomination kind")
        add("shared:" + name, dict(catalogue="sv-candidates-v1", original=target), ends)

    # The website's 1,271-junction denominator collapses calls at the same loci.
    # Preserve its rows independently of our stricter retained-side grouping.
    by_loci = defaultdict(list)
    for number, row in enumerate(fusion_rows, 1):
        loci = ("%s:%s" % (row["chrD"], row["posD"]),
                "%s:%s" % (row["chrA"], row["posA"]))
        name = "fusion-call:%04d" % number
        add(name, dict(catalogue="fusion-survey", original=row), fusion_breakends(row))
        by_loci[loci].append(name)
    survey = []
    for row in summary_rows:
        names = by_loci.get((row["locusA"], row["locusB"]), [])
        if not names:
            raise ValueError("Survey junction has no original sequence call: " + row["label"])
        survey.append(dict(original=row, nominations=names))
    accounted = {name for row in survey for name in row["nominations"]}
    expected = {name for name in targets if name.startswith("fusion-call:")}
    if accounted != expected:
        raise ValueError("Fusion survey and sequence-call inventory disagree")
    return targets, groups, survey


def classify_source(file):
    """Scope from explicit claims and losslessly retained path hints.

    A path hint admits a possible RNA product for inspection, but never resolves
    its biological identity. Blood and unspecified specimens stay separate.
    """
    key = file["key"]
    claims = file["claims"]
    assays = {c["assay"] for c in claims if c.get("assay")}
    published = [c for c in claims if c.get("basis") == "published"]
    tissues = {c["tissue"] for c in published if c.get("tissue")}
    labels = " ".join(c.get("label", "") for c in published).lower()
    hint = bool(re.search(r"(?:^|/)(?:RNA|rna-seq|ONT|pacbio|ucsf|T[123]|hudson_lab)(?:/|$)", key))
    if assays and not assays & RNA_ASSAYS:
        return dict(cohort="excluded", reason="non_RNA_assay")
    if not assays and re.search(r"/(?:wes|wgs|dna)/", key, re.IGNORECASE):
        return dict(cohort="excluded", reason="DNA_path_hint")
    if not assays & RNA_ASSAYS and not hint:
        return dict(cohort="excluded", reason="non_RNA_assay" if assays else "unresolved_assay")
    if re.search(r"(?:vdj[_/]|TCR|BCR|/cider/)", key):
        return dict(cohort="excluded", reason="targeted_immune_repertoire_product")
    if "blood" in tissues or "blood" in labels or re.search(r"(?:/blood/|^hudson_lab/PBMC_)", key):
        cohort = "blood_control"
    elif tissues & {"tumor", "organoid"} or "tumor" in labels or re.search(
            r"(?:^ONT/|^pacbio/|^ucsf/|^T[123]/|/tumor/|^rna-seq/|/RNA/)", key):
        cohort = "tumor_candidate"
    else:
        cohort = "unresolved_specimen"
    if "blood" in tissues and tissues & {"tumor", "organoid"}:
        cohort = "unresolved_specimen"
    return dict(cohort=cohort, reason="published_RNA_assay" if assays & RNA_ASSAYS else "RNA_path_hint",
                assay_claims=sorted(assays), tissue_claims=sorted(tissues),
                timepoint_claims=sorted({c["timepoint"] for c in published if c.get("timepoint")}),
                library_claims=sorted({c["library"] for c in published if c.get("library")}))


def freeze(output, snapshot, cache, metadata_cache):
    """Freeze an existing osteosarc snapshot and explicitly fetch survey tables."""
    import osteosarc

    output = Path(output)
    dataset = osteosarc.Dataset.open(snapshot, cache=cache, offline=True)
    fetcher = osteosarc.Cache(metadata_cache)
    catalogue = osteosarc.load_sv_candidates()
    receipts, tables = {}, {}
    for name in FUSION_TABLES:
        receipt = fetcher.fetch(FUSION_URL + name + ".tsv", max_bytes=10_000_000)
        receipts[name] = receipt.to_dict()
        raw = fetcher.path(receipt).read_bytes()
        tables[name] = list(csv.DictReader(io.StringIO(raw.decode()), delimiter="\t"))
    targets, groups, survey = nominations(catalogue, tables["fusion_sequences"], tables["scan_summary"])
    sources = {}
    for file in dataset.files:
        if file.kind != "alignment":
            continue
        row = asdict(file)
        row.update(selection=classify_source(row), evidence_overlaps=file.evidence_overlaps)
        sources[file.id] = row
    manifest = dict(schema_version=1, assembly="GRCh38", osteosarc_version=osteosarc.__version__,
                    snapshot=dataset.manifest, catalogue_sources=catalogue["sources"],
                    catalogue_identity=identity(catalogue), fusion_receipts=receipts,
                    nominations=targets, geometries=groups, fusion_survey=survey,
                    sources=sources, support_aggregation="within_product_only",
                    curation_diagnostics=list(dataset.corrections),
                    scope="All retained nominations; tumor RNA first, blood controls separately.")
    write_json(output / "inventory.json.gz", manifest)
    write_json(output / "inventory-pin.json", dict(sha256=digest(output / "inventory.json.gz")))
    return manifest


def load_inventory(directory):
    directory = Path(directory)
    pin = read_json(directory / "inventory-pin.json")
    if digest(directory / "inventory.json.gz") != pin["sha256"]:
        raise ValueError("Inventory checksum mismatch")
    return read_json(directory / "inventory.json.gz")


def summary(manifest):
    return dict(nominations=len(manifest["nominations"]), geometries=len(manifest["geometries"]),
                unresolved_nominations=sum(t["geometry_id"] is None for t in manifest["nominations"].values()),
                fusion_survey_junctions=len(manifest["fusion_survey"]),
                sources=dict(Counter(s["selection"]["cohort"] for s in manifest["sources"].values())))
