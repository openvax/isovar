"""Pin Ensembl 87 reference contexts for selected RNA-supported examples.

Uses the existing local PyEnsembl release cache. The PDFs read the resulting
small TSV and do not download or query a live reference during regeneration.
"""

import csv
import gzip
import hashlib
from pathlib import Path

from pyensembl import EnsemblRelease


HERE = Path(__file__).parent
OUTPUT = HERE / "example_reference_contexts.tsv"
EXAMPLES = (
    ("HCC1395", "COSM32759"),
    ("HCC1395", "COSM32783"),
    ("COLO829", "chr9:21971154:CCG>C"),
    ("COLO829", "chr2:105361307:G>GCCCC"),
    ("COLO829", "chr1:156265990:G>A"),
)
FIELDS = ("sample", "call_id", "transcript_id", "protein_start_1based",
          "reference_context", "reference_protein_sha256", "ensembl_release", "placement_method")


def build(output=OUTPUT):
    release = EnsemblRelease(87)
    products = {}
    selected = list(EXAMPLES)
    calls_by_id = {}
    for sample, prefix in (("HCC1395", "hcc1395_direct_rna"),
                           ("COLO829", "colo829")):
        with gzip.open(HERE / f"{prefix}_atlas_proteins.tsv.gz", "rt") as source:
            for row in csv.DictReader(source, delimiter="\t"):
                if row["rank"] == "1":
                    products[sample, row["call_id"]] = row
        with gzip.open(HERE / f"{prefix}_atlas_calls.tsv.gz", "rt") as source:
            calls = sorted(csv.DictReader(source, delimiter='\t'),
                           key=lambda r: (-int(r['alt_fragments']), r['call_id']))
        for category in ('SNV', 'indel'):
            selected.extend((sample, r['call_id']) for r in
                            [r for r in calls if r['class'] == category and (sample, r['call_id']) in products][:3])
        calls_by_id.update({(sample, r['call_id']): r for r in calls})
    rows = []
    for sample, call_id in dict.fromkeys(selected):
        product = products[sample, call_id]
        mutant = product["protein_sequence"]
        start = int(product["mutation_start_idx"])
        end = int(product["mutation_end_idx"])
        for transcript_id in product["compatible_transcript_ids"].split(";"):
            reference = release.transcript_by_id(transcript_id).protein_sequence or ""
            offset = reference.find(mutant[:start])
            if offset < 0 or reference.find(mutant[:start], offset + 1) >= 0:
                continue
            context = reference[offset:offset + len(mutant)]
            if end == start + 1 and context[end:] != mutant[end:]:
                continue
            rows.append(dict(sample=sample, call_id=call_id,
                             transcript_id=transcript_id,
                             protein_start_1based=offset + 1,
                             reference_context=context,
                             reference_protein_sha256=hashlib.sha256(
                                 reference.encode("ascii")).hexdigest(),
                             ensembl_release=87, placement_method='unique_shared_protein_prefix'))
            break
        else:
            call = calls_by_id[sample, call_id]
            for transcript_id in product['compatible_transcript_ids'].split(';'):
                transcript = release.transcript_by_id(transcript_id)
                if not transcript.protein_sequence or not transcript.start_codon_spliced_offsets:
                    continue
                offset = ((transcript.spliced_offset(int(call['pos'])) -
                           min(transcript.start_codon_spliced_offsets)) // 3 - start)
                if offset < 0:
                    continue
                reference = transcript.protein_sequence
                context = reference[offset:offset + len(mutant)]
                if len(context) != len(mutant):
                    continue
                rows.append(dict(sample=sample, call_id=call_id, transcript_id=transcript_id,
                                 protein_start_1based=offset + 1, reference_context=context,
                                 reference_protein_sha256=hashlib.sha256(reference.encode()).hexdigest(),
                                 ensembl_release=87, placement_method='genomic_codon_anchor_RNA_flanks_differ'))
                break
            else:
                raise ValueError(f"No reference context for {sample} {call_id}")
    with output.open("w", newline="") as target:
        writer = csv.DictWriter(target, FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)


if __name__ == "__main__":
    print("Pinned reference contexts:", build())
