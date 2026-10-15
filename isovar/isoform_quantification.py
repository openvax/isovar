"""Salmon abundance estimates kept separate from original RNA witnesses."""

from collections import defaultdict
import csv
from math import isfinite
from pathlib import Path
import subprocess

import pysam

from .audit_bundle import digest, publish_json


def write_isoform_fastas(result, transcript_path, protein_path):
    """Export competitors and complete encoded proteins, retaining provenance IDs.

    Reference competitors can contain uncovered bases and are labelled as such
    in the JSON result. Partial observed sequences have their own result table
    and are not treated as full transcripts during abundance estimation.
    """
    with Path(transcript_path).open("w") as out:
        for row in result["isoforms"]:
            out.write(">%s %s\n%s\n" % (row["id"], row["origin"], row["sequence"]))
    with Path(protein_path).open("w") as out:
        for row in result["isoforms"] + result["partial_sequences"]:
            for i, coding in enumerate(row["coding"]):
                if coding["status"] == "complete_CDS" and all(row["coverage"][coding["query_start"]:coding["query_end"]]):
                    out.write(">%s.protein%d %s\n%s\n" % (row["id"], i + 1, row["origin"], coding["amino_acids"]))


def _indistinguishable(result):
    """Group identical sequences or identical observed fragment compatibility."""
    rows = result["isoforms"]
    parents = {r["id"]: r["id"] for r in rows}

    def root(key):
        while parents[key] != key:
            key = parents[key]
        return key

    seen = {}
    for row in rows:
        cid = row["id"]
        signature = tuple(i for i, c in enumerate(result["support"]["equivalence_classes"])
                          if cid in c["isoform_ids"])
        keys = [("sequence", row["sequence"])]
        if signature:
            keys.append(("evidence", signature))
        for key in keys:
            if key in seen:
                parents[root(cid)] = root(seen[key])
            else:
                seen[key] = cid
    groups = defaultdict(list)
    for cid in parents:
        groups[root(cid)].append(cid)
    return sorted(sorted(ids) for ids in groups.values())


def add_isoform_abundances(result, quant_sf, *, provenance):
    """Attach supplied Salmon estimates in place, without copying read evidence.

    The caller must establish that the quantification used these exact candidate
    sequences among the complete declared competitors. ``quantify_antigen_isoforms``
    builds that input and records the commands/hashes automatically.
    """
    if result.get("schema") != "isovar.antigen_isoforms.v1" or not provenance:
        raise ValueError("A supported isoform result and quantification provenance are required")
    estimates = {}
    with Path(quant_sf).open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            name = row["Name"]
            if name in estimates:
                raise ValueError("Repeated transcript in Salmon estimates: " + name)
            values = {key: float(row[key]) for key in ("Length", "EffectiveLength", "TPM", "NumReads")}
            if any(not isfinite(v) or v < 0 for v in values.values()) or values["Length"] <= 0:
                raise ValueError("Invalid Salmon estimate for " + name)
            estimates[name] = values
    rows = result["isoforms"]
    for row in rows:
        cid = row["id"]
        if cid not in estimates or estimates[cid]["Length"] != len(row["sequence"]):
            raise ValueError("Salmon input does not match reconstructed transcript " + cid)
    groups = _indistinguishable(result)
    genes = []
    for gene in sorted({gene for row in rows for gene in row["gene_ids"]}):
        members = [row for row in rows if gene in row["gene_ids"]]
        ids = {row["id"] for row in members}
        total = sum(estimates[cid]["TPM"] for cid in ids)
        ranked = sorted((dict(isoform_ids=[cid for cid in group if cid in ids],
                             tpm=sum(estimates[cid]["TPM"] for cid in group if cid in ids),
                             also_compatible_with_other_genes=any(cid not in ids for cid in group))
                         for group in groups if ids.intersection(group)), key=lambda g: (-g["tpm"], g["isoform_ids"]))
        for group in ranked:
            group["fraction"] = group["tpm"] / total if total else None
        top = ranked[0] if ranked and total else None
        supported = top is not None and any(row["unique_fragments"] for row in members
                                            if row["id"] in top["isoform_ids"])
        unresolved_paths = any(gene in row["gene_ids"] and not row["compatible_models"]
                               for row in result["partial_sequences"])
        unavailable = any(row["gene_id"] == gene for row in
                          result.get("target_selection", {}).get("unavailable_transcripts", []))
        dominant = (top["isoform_ids"][0] if supported and len(top["isoform_ids"]) == 1
                    and not top["also_compatible_with_other_genes"]
                    and all(row["gene_ids"] == [gene] for row in members if row["id"] in top["isoform_ids"])
                    and (len(ranked) == 1 or top["tpm"] > ranked[1]["tpm"])
                    and not result["discovery"]["truncated"] and not unresolved_paths and not unavailable else None)
        genes.append(dict(gene_id=gene, tpm=total, groups=ranked, dominant_isoform_id=dominant,
                          dominance="estimated" if dominant else "unresolved" if total else "not_observed",
                          isoform_fractions={cid: estimates[cid]["TPM"] / total if total else None
                                             for cid in sorted(ids)}))
    proteins = defaultdict(lambda: dict(isoform_ids=[], tpm=0.0, estimated_fragments=0.0,
                                       observed_coding_tpm=0.0, uncovered_reference_coding_tpm=0.0))
    ambiguous_coding = []
    for row in rows:
        products = {c["amino_acids"] for c in row["coding"] if c["status"] == "complete_CDS"}
        if len(products) != 1 or any(c["status"] != "complete_CDS" for c in row["coding"]):
            ambiguous_coding.append(row["id"])
            continue
        group = proteins[next(iter(products))]
        group["isoform_ids"].append(row["id"])
        group["tpm"] += estimates[row["id"]]["TPM"]
        group["estimated_fragments"] += estimates[row["id"]]["NumReads"]
        observed = all(all(row["coverage"][c["query_start"]:c["query_end"]]) for c in row["coding"])
        group["observed_coding_tpm" if observed else "uncovered_reference_coding_tpm"] += estimates[row["id"]]["TPM"]
    result["quantification"] = dict(
        method="Salmon", provenance=provenance, quant_sf_sha256=digest(quant_sf),
        normalization_scope="supplied_competitive_transcriptome", estimates={r["id"]: estimates[r["id"]] for r in rows},
        background_tpm=sum(value["TPM"] for cid, value in estimates.items() if cid not in {r["id"] for r in rows}),
        indistinguishable_groups=groups, genes=genes,
        proteins=[dict(amino_acids=protein, **row) for protein, row in sorted(proteins.items())],
        unresolved_coding_isoforms=ambiguous_coding,
        interpretation="Model-based RNA abundance; neither protein abundance nor antigen presentation")
    return result


def _augmented_transcriptome(result, source, destination):
    originals = {t["reference"]["transcript_id"]: t["reference"]["sequence"] for t in result["targets"]}
    found, seen = set(), set()
    with pysam.FastxFile(str(source)) as records, Path(destination).open("w") as out:
        for record in records:
            name = record.name.split("|", 1)[0]  # GENCODE FASTA identity field.
            if name in seen or name.startswith("isovar_isoform_"):
                raise ValueError("Duplicate or reserved transcript ID: " + name)
            seen.add(name)
            if name in originals:
                if record.sequence.upper() != originals[name]:
                    raise ValueError("Background transcript sequence differs from selected annotation: " + name)
                found.add(name)
            else:
                out.write(">%s\n%s\n" % (name, record.sequence.upper()))
        if found != originals.keys():
            raise ValueError("Background transcriptome lacks selected annotation IDs: " + ", ".join(sorted(originals.keys() - found)))
        for row in result["isoforms"]:
            out.write(">%s\n%s\n" % (row["id"], row["sequence"]))


def quantify_antigen_isoforms(result, *, transcriptome, reads1, output_dir, reads2=None, salmon="salmon", threads=1):
    """Run installed Salmon on short reads and an augmented full transcriptome.

    The output directory must not exist. No software is installed and no data
    are downloaded. The supplied background must retain all competing genes;
    only exact target reference entries are replaced with the reconstructed set.
    """
    if type(threads) is not int or threads < 1:
        raise ValueError("threads must be a positive integer")
    if result.get("schema") != "isovar.antigen_isoforms.v1":
        raise ValueError("Unsupported isoform result schema")
    output = Path(output_dir).expanduser().resolve()
    inputs = {"transcriptome": Path(transcriptome).expanduser().resolve(),
              "reads1": Path(reads1).expanduser().resolve()}
    if reads2:
        inputs["reads2"] = Path(reads2).expanduser().resolve()
    for path in inputs.values():
        if not path.is_file():
            raise ValueError("Missing quantification input: " + str(path))
    try:
        version = subprocess.run([salmon, "--version"], check=True, capture_output=True, text=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError) as error:
        raise ValueError("A working Salmon installation is required: %s" % error) from error
    output.mkdir(parents=True, exist_ok=False)
    fasta = output / "transcripts.fa"
    _augmented_transcriptome(result, inputs["transcriptome"], fasta)
    index, quant = output / "index", output / "quant"
    commands = [[salmon, "index", "-t", str(fasta), "-i", str(index), "--keepDuplicates", "--no-clip", "-p", str(threads)],
                [salmon, "quant", "-i", str(index), "-l", "A", "-p", str(threads), "-o", str(quant)]]
    commands[1] += (["-1", str(inputs["reads1"]), "-2", str(inputs["reads2"])] if reads2 else
                    ["-r", str(inputs["reads1"])])
    provenance = dict(version=version, commands=commands, transcript_fasta_sha256=digest(fasta),
                      inputs={key: dict(path=str(path), sha256=digest(path)) for key, path in inputs.items()})
    publish_json(output / "inputs.json", provenance)
    for i, command in enumerate(commands):
        with (output / ("index.log" if i == 0 else "quant.log")).open("w") as log:
            try:
                subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
            except (OSError, subprocess.CalledProcessError) as error:
                raise ValueError("Salmon failed; see %s" % log.name) from error
    add_isoform_abundances(result, quant / "quant.sf", provenance=provenance)
    publish_json(output / "isoforms.json", result)
    return result
