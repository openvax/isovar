"""Independent, original-annotation protein expectations for the vaccine cohort.

No Isovar or Varcode consequence code participates in reference selection,
genomic-to-cDNA mapping, editing, or expected translation.
"""

import argparse
from collections import defaultdict
import gzip
from itertools import product
import json
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
from tests.data.osteosarc.expansion.inventory import digest, write_json  # noqa: E402
from tests.data.osteosarc.protein_references import fasta_records  # noqa: E402
from tests.osteosarc_protein_helpers import reverse_complement, transcript_offset  # noqa: E402


# NCBI genetic-code tables 1 and 2, in the published TCAG codon order.
# https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
CODONS = ["".join(b) for b in product("TCAG", repeat=3)]
TABLES = {
    1: dict(zip(CODONS, "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")),
    2: dict(zip(CODONS, "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG")),
}
STARTS = {1: {"TTG", "CTG", "ATG"}, 2: {"ATT", "ATC", "ATA", "ATG", "GTG"}}


def translate(sequence, table=1, annotated_start=False):
    amino_acids = []
    for offset in range(0, len(sequence) - 2, 3):
        codon = sequence[offset:offset + 3]
        aa = TABLES[table].get(codon, "X")
        if offset == 0 and annotated_start and codon in STARTS[table]:
            aa = "M"
        if aa == "*":
            return "".join(amino_acids), True
        amino_acids.append(aa)
    return "".join(amino_acids), False


def minimal_edit(record):
    """Remove identical sequence only; do not repeat-normalize distinct loci."""
    pos, ref, alt = int(record["pos"]), record["ref"], record["alt"]
    while ref and alt and ref[0] == alt[0]:
        pos, ref, alt = pos + 1, ref[1:], alt[1:]
    while ref and alt and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    if not ref and not alt:
        raise ValueError("Identical reference and alternate are not a variant")
    return pos, ref, alt


def apply_variant(record, model):
    """Map and edit an exact exonic allele, including reverse-strand insertions."""
    if record["chrom"].removeprefix("chr").replace("MT", "M") != model["contig"].replace("MT", "M"):
        raise ValueError("Variant/reference contig mismatch")
    pos, ref, alt = minimal_edit(record)
    if ref:
        offsets = sorted(transcript_offset(model["exons"], model["strand"], p)
                         for p in range(pos, pos + len(ref)))
        if offsets != list(range(offsets[0], offsets[0] + len(ref))):
            raise ValueError("Variant is not contiguous in this transcript")
        offset = offsets[0]
    else:
        left, right = [transcript_offset(model["exons"], model["strand"], p) for p in (pos - 1, pos)]
        if abs(left - right) != 1:
            raise ValueError("Insertion anchors are not adjacent in this transcript")
        offset = max(left, right)
    if model["strand"] == "-":
        ref, alt = reverse_complement(ref), reverse_complement(alt)
    if model["cdna"][offset:offset + len(ref)] != ref:
        raise ValueError("Reference allele does not match original transcript cDNA")
    if offset < model["cds_start"] or offset >= model["cds_end"]:
        raise ValueError("Variant is outside the annotated coding region")
    mutant_cdna = model["cdna"][:offset] + alt + model["cdna"][offset + len(ref):]
    mutant_protein, stop = translate(mutant_cdna[model["cds_start"]:], model["genetic_code"], True)
    return dict(model, variant_offset=offset, oriented_ref=ref, oriented_alt=alt,
                mutant_cdna=mutant_cdna, mutant_protein=mutant_protein,
                mutant_has_stop=stop, frameshift=(len(ref) - len(alt)) % 3 != 0)


def check_translation(translation, expected, protein_length):
    """Check RNA translation and its independent nominated-edit expectation."""
    orf = translation.variant_orf
    prefix = orf.variant_cdna_interval_start
    start = expected["variant_offset"] - prefix
    if start < 0:
        raise ValueError("RNA context extends before this annotated transcript")
    frame = (expected["cds_start"] - start if start < expected["cds_start"]
             else (expected["cds_start"] - start) % 3)
    assert orf.offset_to_first_complete_codon == frame, "reading-frame mismatch"
    assert len(orf.reference_cdna_sequence_before_variant) == prefix, "reference-prefix offset mismatch"
    assert orf.cdna_sequence[prefix:orf.variant_cdna_interval_end] == expected["oriented_alt"], "alternate allele mismatch"
    assert translation.frameshift == expected["frameshift"], "frameshift flag mismatch"
    annotated_start = start + frame == expected["cds_start"]
    actual, actual_stop = translate(orf.cdna_sequence[frame:], expected["genetic_code"], annotated_start)
    cdna = expected["mutant_cdna"][start:start + len(orf.cdna_sequence)]
    reference, reference_stop = translate(cdna[frame:], expected["genetic_code"], annotated_start)
    if protein_length and len(actual) > protein_length:
        actual, actual_stop = actual[:protein_length], False
    if protein_length and len(reference) > protein_length:
        reference, reference_stop = reference[:protein_length], False
    assert (translation.amino_acids, translation.ends_with_stop_codon) == (actual, actual_stop), "RNA translation/stop mismatch"
    protein_start = (start + frame - expected["cds_start"]) // 3
    mutant_start = (expected["variant_offset"] - expected["cds_start"]) // 3 - protein_start
    mutant_end = (len(actual) if expected["frameshift"] else
                  (expected["variant_offset"] + len(expected["oriented_alt"]) - expected["cds_start"] + 2) // 3 - protein_start)
    assert translation.mutation_start_idx == mutant_start, "mutation start mismatch"
    assert translation.mutation_end_idx == (min(mutant_end, protein_length) if protein_length else mutant_end), "mutation end mismatch"
    contains = (0 < mutant_start < len(actual) if mutant_start == mutant_end else len(actual) > mutant_start)
    assert translation.contains_mutation == contains, "mutation containment mismatch"
    return dict(transcript_version=expected["transcript_version"], protein_version=expected["protein_version"],
                protein_start_1based=protein_start + 1, expected_amino_acids=reference,
                expected_stop=reference_stop, matches_expected=(actual, actual_stop) == (reference, reference_stop),
                frame=frame, genetic_code=expected["genetic_code"])


def load_reference(directory):
    """Verify every offline asset before using any derived expectation."""
    directory = Path(directory)
    manifest = json.loads((directory / "manifest.json").read_text())
    for name, checksum in manifest["files"].items():
        if digest(directory / name) != checksum:
            raise ValueError(f"Reference asset checksum mismatch: {name}")
    with gzip.open(directory / "models.json.gz", "rt") as handle:
        models = json.load(handle)
    return manifest, models


def reference_genome(directory, cache):
    """Use real PyEnsembl parsing with an explicit, distinct subset identity."""
    from pyensembl import Genome

    directory = Path(directory)
    manifest, _ = load_reference(directory)
    genome = Genome(
        reference_name=manifest["dataset_identity"], annotation_name="osteosarc-vaccine-cohort",
        annotation_version=manifest["ensembl_release"], gtf_path_or_url=str(directory / "reference.gtf.gz"),
        transcript_fasta_paths_or_urls=[str(directory / "reference.cdna.fa.gz")],
        protein_fasta_paths_or_urls=[str(directory / "reference.pep.fa.gz")],
        copy_local_files_to_cache=True, cache_directory_path=str(cache))
    genome.index()
    return genome


def gzip_bytes(path, value):
    with Path(path).open("xb") as handle:
        with gzip.GzipFile(filename="", mode="wb", fileobj=handle, mtime=0) as archive:
            archive.write(value)


def raw_features(path):
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                yield line, None, None
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9:
                raise ValueError("Malformed original GTF row")
            yield line, fields, dict(re.findall(r'(\S+) "([^"]*)";', fields[8]))


def build_reference(source, output, variants, assembly="GRCh38", release=87):
    source, output = Path(source), Path(output)
    output.mkdir(parents=True, exist_ok=False)
    files = {
        "gtf": source / f"Homo_sapiens.{assembly}.{release}.gtf.gz",
        "cdna": source / f"Homo_sapiens.{assembly}.cdna.all.fa.gz",
        "pep": source / f"Homo_sapiens.{assembly}.pep.all.fa.gz",
    }
    # Older Ensembl archives include the release number in FASTA filenames.
    for kind in ("cdna", "pep"):
        if not files[kind].exists():
            files[kind] = source / f"Homo_sapiens.{assembly}.{release}.{kind}.all.fa.gz"
    by_contig = defaultdict(list)
    for record in variants:
        by_contig[record["chrom"].removeprefix("chr").replace("MT", "M")].append(record)
    tids, genes = set(), set()
    for _, fields, attrs in raw_features(files["gtf"]):
        if fields is None or fields[2] != "transcript":
            continue
        records = by_contig[fields[0].replace("MT", "M")]
        if any(int(fields[3]) <= r["pos"] <= int(fields[4]) for r in records):
            tids.add(attrs["transcript_id"])
            genes.add(attrs["gene_id"])
    features, lines = defaultdict(lambda: defaultdict(list)), []
    for line, fields, attrs in raw_features(files["gtf"]):
        if (fields is None or attrs.get("transcript_id") in tids
                or (fields[2] == "gene" and attrs["gene_id"] in genes)):
            lines.append(line)
            if fields is not None and "transcript_id" in attrs:
                features[attrs["transcript_id"]][fields[2]].append(fields)
    cdnas = {h[1:].split()[0].split(".")[0]: (h, s) for h, s in fasta_records(files["cdna"])
             if h[1:].split()[0].split(".")[0] in tids}
    peptides = {re.search(r"transcript:(ENST\d+)", h)[1]: (h, s)
                for h, s in fasta_records(files["pep"])
                if re.search(r"transcript:(ENST\d+)", h)[1] in tids}
    models, unavailable = {}, {}
    for tid in sorted(tids):
        rows = features[tid]
        if tid not in cdnas or tid not in peptides or not rows["start_codon"] or not rows["CDS"]:
            unavailable[tid] = "incomplete_or_noncoding_reference_model"
            continue
        strand = rows["transcript"][0][6]
        exons = [(int(r[3]), int(r[4])) for r in rows["exon"]]
        starts = [transcript_offset(exons, strand, p) for r in rows["start_codon"]
                  for p in range(int(r[3]), int(r[4]) + 1)]
        cds_start = min(starts)
        cds_end = max(transcript_offset(exons, strand, p) for r in rows["CDS"]
                      for p in (int(r[3]), int(r[4]))) + 1
        header, cdna = cdnas[tid]
        protein_header, protein = peptides[tid]
        if sum(e - s + 1 for s, e in exons) != len(cdna):
            unavailable[tid] = "original_cdna_exon_length_mismatch"
            continue
        table = 2 if rows["transcript"][0][0] == "MT" else 1
        translated, has_stop = translate(cdna[cds_start:], table, True)
        if translated != protein:
            unavailable[tid] = "original_cdna_peptide_translation_mismatch"
            continue
        models[tid] = dict(
            transcript_id=tid, transcript_version=header[1:].split()[0],
            protein_version=protein_header[1:].split()[0], contig=rows["transcript"][0][0],
            exons=exons, strand=strand, cds_start=cds_start, cds_end=cds_end,
            cdna=cdna, protein=protein, genetic_code=table, reference_has_stop=has_stop,
        )
    variant_models, variant_gaps = {}, {}
    for record in variants:
        matches, gaps = [], {}
        for tid, model in models.items():
            if model["contig"].replace("MT", "M") != record["chrom"].removeprefix("chr").replace("MT", "M"):
                continue
            if not min(s for s, _ in model["exons"]) <= record["pos"] <= max(e for _, e in model["exons"]):
                continue
            try:
                apply_variant(record, model)
                matches.append(tid)
            except ValueError as error:
                gaps[tid] = str(error)
        variant_models[record["variant_id"]] = matches
        variant_gaps[record["variant_id"]] = gaps
    assets = {"reference.gtf.gz": "".join(lines).encode()}
    for name, records in (("reference.cdna.fa.gz", cdnas), ("reference.pep.fa.gz", peptides)):
        assets[name] = "".join(h + "\n" + s + "\n" for h, s in records.values()).encode()
    for name, content in assets.items():
        gzip_bytes(output / name, content)
    gzip_bytes(output / "models.json.gz", json.dumps(models, sort_keys=True).encode())
    manifest = dict(
        assembly=assembly, ensembl_release=release,
        dataset_identity=f"{assembly}-osteosarc-vaccine-cohort-ensembl{release}",
        sources={kind: dict(filename=path.name, sha256=digest(path),
                            url=f"https://ftp.ensembl.org/pub/release-{release}/"
                                + ("gtf/homo_sapiens/" if kind == "gtf" else f"fasta/homo_sapiens/{kind}/") + path.name)
                 for kind, path in files.items()},
        files={name: digest(output / name) for name in [*assets, "models.json.gz"]},
        transcript_count=len(models), variant_transcripts=variant_models,
        unavailable_transcripts=unavailable, variant_transcript_gaps=variant_gaps,
        scope="all independently validated complete coding models at the pinned variant loci; explicit model exclusions",
    )
    write_json(output / "manifest.json", manifest)
    print("REFERENCE", assembly, "models", len(models), "variants", sum(bool(x) for x in variant_models.values()),
          "of", len(variants), "unavailable models", len(unavailable), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("inventory", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--assembly", default="GRCh38")
    parser.add_argument("--release", type=int, default=87)
    args = parser.parse_args()
    build_reference(args.source, args.output, json.loads(args.inventory.read_text())["variants"],
                    args.assembly, args.release)
