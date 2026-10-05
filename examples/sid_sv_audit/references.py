"""Pin all transcripts of genes intersecting either breakend neighbourhood."""

from pathlib import Path

from isovar.fusion import FusionReference
from .inventory import digest, load_inventory, read_json, write_json


def build(directory, release=115, window=2000):
    from pyensembl import EnsemblRelease

    directory = Path(directory)
    manifest = load_inventory(directory)
    genome = EnsemblRelease(release)
    models, assignments, excluded, cds_unavailable = {}, {}, {}, {}
    for name, group in manifest["geometries"].items():
        gene_ids = set()
        for end in group["breakends"]:
            # PyEnsembl uses one-based inclusive loci; the input is interbase.
            contig, pos = end["contig"].removeprefix("chr"), end["position"]
            gene_ids.update(genome.gene_ids_at_locus(contig, max(1, pos - window + 1), pos + window))
        transcript_ids = {tid for gid in gene_ids for tid in genome.transcript_ids_of_gene_id(gid)}
        for tid in sorted(transcript_ids - models.keys() - excluded.keys()):
            t = genome.transcript_by_id(tid)
            sequence = t.sequence
            if not sequence or set(sequence) - set("ACGT"):
                excluded[tid] = "missing_or_non_ACGT_transcript_sequence"
                continue
            cds_start = min(t.start_codon_spliced_offsets) if t.complete else None
            cds_end = max(t.stop_codon_spliced_offsets) + 1 if t.complete else None
            model = dict(transcript_id=tid, reference_name="GRCh38", annotation="Ensembl %d" % release,
                         contig="chr" + t.contig, strand=t.strand,
                         exons=sorted((a - 1, b) for a, b in t.exon_intervals),
                         sequence=sequence, cds_start=cds_start, cds_end=cds_end)
            # Validate the same contract the reconstruction engine will use.
            try:
                FusionReference(**model)
            except ValueError as error:
                # Ensembl also supplies noncanonical or recoded CDSs. Preserve
                # their sequence/exons for splice and UTR comparisons without
                # inventing a standard-code frame. Retain the original CDS.
                if cds_start is None:
                    raise
                unframed = dict(model, cds_start=None, cds_end=None)
                FusionReference(**unframed)
                cds_unavailable[tid] = dict(reason=str(error), cds_start=cds_start, cds_end=cds_end,
                                            sequence=sequence[cds_start:cds_end])
                model = unframed
            models[tid] = model
        assignments[name] = dict(transcripts=sorted(transcript_ids & models.keys()),
                                 excluded_transcripts=sorted(transcript_ids & excluded.keys()),
                                 gene_ids=sorted(gene_ids))
    paths = [genome.gtf_path, *genome.transcript_fasta_paths]
    result = dict(annotation="Ensembl %d" % release, window=window,
                  inventory_sha256=digest(directory / "inventory.json.gz"),
                  source_files=[dict(name=Path(p).name, sha256=digest(p)) for p in paths],
                  models=models, assignments=assignments, excluded=excluded,
                  unsupported_cds=cds_unavailable)
    write_json(directory / "references.json.gz", result)
    write_json(directory / "references-pin.json", dict(sha256=digest(directory / "references.json.gz")))
    return result


def load(directory):
    directory = Path(directory)
    pin = read_json(directory / "references-pin.json")
    if digest(directory / "references.json.gz") != pin["sha256"]:
        raise ValueError("Reference checksum mismatch")
    result = read_json(directory / "references.json.gz")
    if result["inventory_sha256"] != digest(directory / "inventory.json.gz"):
        raise ValueError("References belong to a different inventory")
    return result
