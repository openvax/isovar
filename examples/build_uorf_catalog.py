"""Rebuild the offline uORF reference from checksum-pinned public source data.

Run as ``python -m examples.build_uorf_catalog --help`` from the checkout.
Only building needs openpyxl; preparing the author's Rdata also needs Rscript.
The installed catalogue and query API need neither, and never access the network.
"""

import argparse
from collections import defaultdict
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import subprocess
import tarfile
import urllib.request
import zipfile

from isovar.uorf_catalog import load_uorf_catalog


COMMIT = "b65ae0c09a53e9efe083d806d617581e6139bdde"
AUTHOR = "https://raw.githubusercontent.com/VanHeeschLab/deutsch_kok_et_al_2024/" + COMMIT + "/raw/"
PINS = {
    "ncorf_list.xlsx": "ca2e2f9389052c62f48ec6f81d53ca88c64152edb5940e556813d4a7bf1e6c62",
    "filtered_peptides.tar.bz2": "a97d73b08b4ee7376a064f9fee3480b3d73639ebe0bd238a6dd0099df6a9e30c",
    "peptide_sample_msrun_counts.tar.bz2": "290658de32de2579d17f0b75c3fcc39b4d50c5e6ce1ceea9ce9f411ccc8a98ba",
    "author-hla-mappings.tsv": "7f97f4ea675ac6abc97217c865019d73135358f448533ad8d98c9479fee9fb8f",
    "author-hla-counts.json": "b1e8e79e74e90d205224cc86663ebf9d64ed97cf0e0290cd86117756ee46a4b6",
    "media-4_with_dot_products.csv": "9e2cbb54765f8a9886cdcbf0bdfefa772f4623f54291add2bf3c47c26b17250f",
}
RIBOSOME = {
    "Ji_et_al_2015": "26687005", "Calviello_et_al_2016": "26657557",
    "Raj_et_al_2016": "27232982", "VanHeesch_et_al_2019": "31155234",
    "Martinez_et_al_2020": "31819274", "Chen_et_al_2020": "32139545",
    "Gaertner_et_al_2020": "32744504",
}
# These legacy aggregate annotations establish whole-proteome MS/MS only.
# Mixed HLA/whole-proteome or targeted/shotgun studies remain unclassified.
SHOTGUN_STUDIES = {
    "Vanderperre_et_al_2013", "Kim_et_al_2014", "Calviello_et_al_2016",
    "Raj_et_al_2016", "Zhu_et_al_2018", "Cao_et_al_2020",
}


def digest(path):
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def verify(path):
    if digest(path) != PINS[path.name]:
        raise ValueError("Source checksum mismatch: %s" % path)


class RemoteZip(io.RawIOBase):
    """Read small ZIP members with HTTP ranges, never the 4.8 GB archive."""

    size = 4797315346
    url = "https://ndownloader.figshare.com/files/57979642"

    def __init__(self):
        self.position = 0

    def seekable(self):
        return True

    def tell(self):
        return self.position

    def seek(self, offset, whence=0):
        self.position = offset if whence == 0 else self.position + offset if whence == 1 else self.size + offset
        return self.position

    def read(self, count=-1):
        if count < 0:
            count = self.size - self.position
        if not 0 <= count <= 2_000_000:
            raise ValueError("Refusing a large archive range")
        if not count:
            return b""
        request = urllib.request.Request(self.url, headers={
            "Range": "bytes=%d-%d" % (self.position, self.position + count - 1)})
        with urllib.request.urlopen(request, timeout=120) as response:
            if response.status != 206:
                raise ValueError("Source server does not support byte ranges")
            data = response.read(count)
        self.position += len(data)
        return data


def prepare_sources(directory):
    """Download ~135 MB; extract small tables and stream the count archive."""
    directory.mkdir(parents=True, exist_ok=True)
    for name in ("ncorf_list.xlsx", "filtered_peptides.tar.bz2", "peptide_sample_msrun_counts.tar.bz2"):
        target = directory / name
        if not target.exists():
            with urllib.request.urlopen(AUTHOR + name, timeout=120) as response, target.open("wb") as out:
                while block := response.read(1024 * 1024):
                    out.write(block)
        verify(target)
    mappings = directory / "author-hla-mappings.tsv"
    if not mappings.exists():
        # The source Rdata is ~825 MiB. Load in an isolated R environment,
        # export only the small HLA mapping table, then remove this scratch file.
        rdata = directory / "filtered_peptides.Rdata"
        with tarfile.open(directory / "filtered_peptides.tar.bz2", "r:bz2") as archive:
            with archive.extractfile("filtered_peptides.Rdata") as source, rdata.open("wb") as out:
                while block := source.read(1024 * 1024):
                    out.write(block)
        try:
            script = ('e<-new.env();load(commandArgs(TRUE)[1],envir=e);'
                      'write.table(e$pep_nonc_long,commandArgs(TRUE)[2],sep="\\t",'
                      'quote=FALSE,row.names=FALSE,na="")')
            subprocess.run(["Rscript", "-e", script, str(rdata), str(mappings)], check=True)
        finally:
            rdata.unlink()
    verify(mappings)
    counts_path = directory / "author-hla-counts.json"
    if not counts_path.exists():
        with mappings.open() as stream:
            peptides = {r["peptide"] for r in csv.DictReader(stream, delimiter="\t")}
        counts, datasets, runs = defaultdict(int), defaultdict(set), defaultdict(set)
        with tarfile.open(directory / "peptide_sample_msrun_counts.tar.bz2", "r|bz2") as archive:
            for member in archive:
                if member.name.endswith("peptide_sample_msrun_counts.tsv"):
                    for raw in archive.extractfile(member):
                        peptide, sample, count, dataset, tag, run = raw.decode().rstrip("\n").split("\t")
                        if peptide in peptides:
                            counts[peptide] += int(count)
                            datasets[peptide].add(dataset)
                            runs[peptide].add((sample, run))
        result = {p: dict(spectrum_count=counts[p], datasets=sorted(datasets[p]), run_count=len(runs[p]))
                  for p in sorted(counts)}
        counts_path.write_text(json.dumps(result, sort_keys=True, separators=(",", ":")) + "\n")
    verify(counts_path)
    benchmark = directory / "media-4_with_dot_products.csv"
    if not benchmark.exists():
        with zipfile.ZipFile(RemoteZip()) as archive:
            benchmark.write_bytes(archive.read("input_files/media-4_with_dot_products.csv"))
    verify(benchmark)


def split(value):
    return tuple(x for x in str(value or "").split(";") if x)


def intervals(protein, peptide):
    protein, peptide = protein.replace("I", "L"), peptide.replace("I", "L")
    result, start = [], 0
    while (position := protein.find(peptide, start)) >= 0:
        result.append([position + 1, position + len(peptide)])
        start = position + 1
    return result


def peptide_record(model, peptide, **values):
    locations = intervals(model["protein_sequence"], peptide)
    if not locations:
        raise ValueError("Peptide absent from %s: %s" % (model["orf_id"], peptide))
    return dict(sequence=peptide, protein_intervals=locations,
                sequence_match="exact" if all(model["protein_sequence"][a - 1:b] == peptide
                                              for a, b in locations) else "I_L_equivalent",
                **values)


def build(source_dir, destination):
    import openpyxl

    for name in ("ncorf_list.xlsx", "author-hla-mappings.tsv", "author-hla-counts.json", "media-4_with_dot_products.csv"):
        verify(source_dir / name)
    workbook = openpyxl.load_workbook(source_dir / "ncorf_list.xlsx", read_only=True, data_only=True)
    models, aliases, excluded = {}, {}, []
    for sheet in ("S2. PHASE I Ribo-seq ORFs", "S3. Single-study Ribo-seq ORFs"):
        rows = workbook[sheet].iter_rows(values_only=True)
        header = next(rows)
        for values in rows:
            row = dict(zip(header, values))
            if row.get("orf_biotype") not in ("uORF", "uoORF"):
                continue
            protein = row["orf_sequence"]
            if not protein.endswith("*") or "*" in protein[:-1]:
                raise ValueError("Unexpected unterminated ORF: %s" % row["orf_name"])
            model = dict(orf_id=row["orf_name"], gene_id=row["gene_id"], gene_name=row["gene_name"],
                         transcript_id=row["transcript"], contig=str(row["chrm"]), strand=row["strand"],
                         coding_blocks=[list(x) for x in zip(map(int, split(row["starts"])), map(int, split(row["ends"])))],
                         protein_sequence=protein[:-1], biotype=row["orf_biotype"], genome_build="GRCh38",
                         ribosome_studies=[s for s in RIBOSOME if row[s] == 1], peptides=[])
            models[model["orf_id"]] = model
            aliases[model["orf_id"]] = set(split(row["literature_ORF_names"])) | {model["orf_id"]}
            datasets = sorted(set(split(row["MS_datasets"])))
            method = "shotgun" if datasets and set(datasets) <= SHOTGUN_STUDIES else "ms_unclassified"
            for peptide in sorted(set(split(row["MS_peptides"]))):
                if not intervals(model["protein_sequence"], peptide):
                    excluded.append(dict(orf_id=model["orf_id"], sequence=peptide,
                                         source="gencode_phase1", reason="sequence_not_in_exact_model"))
                    continue
                model["peptides"].append(peptide_record(model, peptide, method=method,
                    source="gencode_phase1", accession=model["orf_id"], datasets=datasets))
    workbook.close()
    counts = json.loads((source_dir / "author-hla-counts.json").read_text())
    with (source_dir / "author-hla-mappings.tsv").open() as stream:
        seen = set()
        for row in csv.DictReader(stream, delimiter="\t"):
            model = models.get(row["orf_name"])
            key = (row["orf_name"], row["peptide_accession"])
            if model is None or key in seen:
                continue
            seen.add(key)
            if model["protein_sequence"].replace("I", "L") != row["protein"].replace("I", "L"):
                excluded.append(dict(orf_id=model["orf_id"], sequence=row["peptide"],
                                     source="author_hla", reason="source_protein_differs_from_exact_model"))
                continue
            count = counts.get(row["peptide"], {})
            model["peptides"].append(peptide_record(model, row["peptide"], method="hla",
                source="author_hla", accession=row["peptide_accession"],
                mapping_count=int(row["n_pep"]), spectrum_count=count.get("spectrum_count"),
                datasets=count.get("datasets", [])))
    # Keep each spectrum and its reviewers together. A compatible sequence in
    # another isoform is explicit, and never receives a uniqueness claim.
    reviews = defaultdict(list)
    with (source_dir / "media-4_with_dot_products.csv").open(encoding="utf-8-sig") as stream:
        for row in csv.DictReader(stream):
            if row["pmid"] != "Control" and row["Rating..1.5."] not in ("", "NA"):
                reviews[row["USI"]].append(row)
    for usi, rows in sorted(reviews.items()):
        first = rows[0]
        scores = sorted((r["curator"], int(r["Rating..1.5."])) for r in rows)
        quality = ("reviewed_support" if all(s >= 4 for _, s in scores) else
                   "reviewed_low_quality" if all(s < 4 for _, s in scores) else "reviewed_mixed")
        for model in models.values():
            if intervals(model["protein_sequence"], first["peptide"]):
                model["peptides"].append(peptide_record(model, first["peptide"],
                    method="hla" if first["HLA"] == "TRUE" else "shotgun",
                    source="microprotein_benchmark", accession=usi, spectrum_ids=[usi],
                    quality=quality, review_ratings=scores, datasets=[usi.split(":")[1]],
                    attribution="source_orf" if first["orf"] in aliases[model["orf_id"]] else "sequence_compatible"))
    for model in models.values():
        model["peptides"].sort(key=lambda p: (p["source"], p["sequence"], p["accession"]))
    metadata = dict(dataset_version="2026-10-08.1", genome_build="GRCh38", coordinates="1-based-inclusive",
        model_scope="3771 stopped uORF/uoORF models from the GENCODE Phase I 2022 reference (minimum 16 aa)",
        evidence_scope="Legacy reported MS, author HLA mapping export, and separately licensed community spectrum reviews; not the complete 2026 Nature curated shotgun tables",
        missing_evidence="not_assessed_or_not_reported; never evidence of absence",
        mutant_translation_evidence="not_assessed",
        benchmark_review_policy="All reviewers score 4 or 5: reviewed_support; all below 4: reviewed_low_quality; otherwise reviewed_mixed",
        ribosome_study_pmids=RIBOSOME,
        sources=[
            dict(id="gencode_phase1", citation="Mudge et al. 2022, Standardized annotation of translated open reading frames",
                 doi="10.1038/s41587-022-01369-0", url=AUTHOR + "ncorf_list.xlsx", sha256=PINS["ncorf_list.xlsx"],
                 license="MIT (author repository)", mapping_scope="No peptide uniqueness counts; legacy MS datasets are ORF-level pooled annotations, not peptide-specific sample attribution"),
            dict(id="author_hla", citation="Deutsch et al. 2026, Expanding the human proteome with microproteins and peptideins",
                 doi="10.1038/s41586-026-10459-x", url=AUTHOR + "filtered_peptides.tar.bz2", sha256=PINS["filtered_peptides.tar.bz2"],
                 counts_url=AUTHOR + "peptide_sample_msrun_counts.tar.bz2", counts_sha256=PINS["peptide_sample_msrun_counts.tar.bz2"],
                 license="MIT (author repository)", mapping_scope="n_pep in author's filtered HLA reference mappings, excluding canonical peptides; not global proteome uniqueness",
                 count_scope="Peptide PSMs across source MS runs, shared across ORF mappings; not RNA abundance or independent patient samples"),
            dict(id="microprotein_benchmark", citation="Wacholder et al. 2026, Community benchmarking and evaluation of human unannotated microprotein detection by mass spectrometry based proteomics",
                 doi="10.1038/s41467-025-68002-x", dataset_doi="10.6084/m9.figshare.30131869.v1",
                 url=RemoteZip.url, member="input_files/media-4_with_dot_products.csv", sha256=PINS["media-4_with_dot_products.csv"],
                 license="CC-BY-4.0 (Figshare dataset)", mapping_scope="Sequence-compatible Phase I models; source_orf only when a published ORF alias also matches; uniqueness unknown",
                 correction="10.1038/s41467-026-73431-3 corrects figure labels, not source spectra")])
    data = dict(schema="isovar.uorf_evidence.v1", metadata=metadata, records=sorted(models.values(), key=lambda m: m["orf_id"]))
    destination.mkdir(parents=True, exist_ok=True)
    target = destination / "catalog.json.gz"
    target.write_bytes(gzip.compress((json.dumps(data, sort_keys=True, separators=(",", ":")) + "\n").encode(), mtime=0))
    catalog = load_uorf_catalog(target)
    stats = dict(records=len(catalog.records), peptides=sum(len(r.peptides) for r in catalog.records),
                 excluded_legacy_assignments=len(excluded),
                 **{method + "_orfs": len(catalog.query(evidence=method)) for method in ("ribosome", "shotgun", "hla", "ms")},
                 reviewed_ms_orfs=len(catalog.query(evidence="ms", require_reviewed=True)))
    manifest = dict(schema="isovar.uorf_manifest.v1", dataset_version=metadata["dataset_version"],
                    catalog_sha256=digest(target), statistics=stats, source_sha256=PINS,
                    excluded_assignments=excluded)
    (destination / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(stats, sort_keys=True))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--prepare-sources", action="store_true", help="Download checksum-pinned inputs; requires Rscript")
    parser.add_argument("--output", type=Path, default=Path("isovar/data/uorf-evidence"))
    args = parser.parse_args()
    if args.prepare_sources:
        prepare_sources(args.source_dir)
    build(args.source_dir, args.output)


if __name__ == "__main__":
    main()
