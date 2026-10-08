"""Query and export the offline uORF evidence reference."""

import argparse
import csv
import json
import re
import sys

from ..uorf_catalog import EVIDENCE_METHODS, load_uorf_catalog


def run(argv=None, prog=None):
    """Export reference records as JSON, TSV or FASTA."""
    parser = argparse.ArgumentParser(description=__doc__, prog=prog)
    parser.add_argument("--catalog", help="Use a compatible local JSON/JSON.gz catalogue")
    parser.add_argument("--gene", help="Exact gene symbol or Ensembl gene ID")
    parser.add_argument("--transcript", help="Exact source transcript ID (unversioned)")
    parser.add_argument("--orf-id", help="Exact catalogue ORF identifier")
    parser.add_argument("--evidence", choices=EVIDENCE_METHODS)
    parser.add_argument("--require-unique", action="store_true", help="Require one protein mapping in the MS source")
    parser.add_argument("--require-reviewed", action="store_true", help="Require convincing reviewed spectra (benchmark scores 4-5)")
    parser.add_argument("--biotype", choices=("uORF", "uoORF"))
    parser.add_argument("--region", help="1-based inclusive interval: chr7:66205483-66205522")
    parser.add_argument("--position", help="Map a 1-based variant position: chr7:66205483 (JSON only)")
    parser.add_argument("--genome-build", help="Required for genomic queries; packaged data use GRCh38")
    parser.add_argument("--format", choices=("json", "tsv", "fasta"), default="json")
    args = parser.parse_args(argv)
    if args.region and args.position:
        parser.error("Choose --region or --position")
    if args.position and args.format != "json":
        parser.error("--position requires --format json")
    try:
        catalog = load_uorf_catalog(args.catalog)
        options = dict(gene=args.gene, transcript_id=args.transcript,
                       evidence=args.evidence, require_unique=args.require_unique,
                       require_reviewed=args.require_reviewed,
                       biotype=args.biotype, genome_build=args.genome_build)
        locus = args.region or args.position
        if locus:
            pattern = r"([^:]+):([0-9]+)-([0-9]+)" if args.region else r"([^:]+):([0-9]+)"
            match = re.fullmatch(pattern, locus)
            if match is None:
                raise ValueError("Invalid locus: %s" % locus)
            start = int(match.group(2))
            options.update(contig=match.group(1), start=start,
                           end=int(match.group(3)) if args.region else start)
        records = catalog.query(**options)
        if args.orf_id:
            catalog.get(args.orf_id)  # A misspelled ID is an error, not empty evidence.
            records = tuple(r for r in records if r.orf_id == args.orf_id)
    except (ValueError, KeyError, OSError) as exc:
        parser.error(str(exc))
    if args.format == "json":
        output = dict(schema="isovar.uorf_query.v1", metadata=catalog.metadata,
                      records=[r.to_dict() for r in records])
        if args.position:
            output["position_annotations"] = [
                r.map_genomic_position(options["start"], genome_build=args.genome_build)
                for r in records]
        json.dump(output, sys.stdout, indent=2)
        sys.stdout.write("\n")
    elif args.format == "fasta":
        for r in records:
            sys.stdout.write(">%s gene=%s transcript=%s assembly=%s biotype=%s\n%s\n" % (
                r.orf_id, r.gene_name, r.transcript_id, r.genome_build, r.biotype, r.protein_sequence))
    else:
        fields = ("orf_id", "gene_id", "gene_name", "transcript_id", "genome_build",
                  "contig", "strand", "biotype", "protein_sequence", "coding_blocks_json",
                  "ribosome_studies_json", "peptides_json")
        writer = csv.DictWriter(sys.stdout, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for r in records:
            row = r.to_dict()
            for field in ("coding_blocks", "ribosome_studies", "peptides"):
                row[field + "_json"] = json.dumps(row.pop(field), separators=(",", ":"))
            writer.writerow(row)


if __name__ == "__main__":
    run()
