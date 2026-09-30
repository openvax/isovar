"""Independent CIGAR evidence and shared inputs for the full Sid panel.

The allele oracle uses original SAM operations and literal target alleles;
it does not invoke Isovar's collection or Varcode's normalization.
"""

from collections import Counter

from tests.data.osteosarc.expansion.references import minimal_edit
from tests.real_rna_helpers import cigar_observation, record_digest


def target_record(name, target):
    return dict(variant_id=name, chrom=target["contig"], pos=target["position"],
                ref=target["ref"], alt=target["alt"])


def direct_observation(read, record):
    """Aligned allele for SNVs, or two flanking anchors for other literal edits.

    Complex replacements and MNVs use the minimal changed interval. N skips
    cannot support deletions. Unanchored indels remain uncallable here; this
    oracle intentionally does not infer repeat-equivalent indels from MD.
    """
    if len(record["ref"]) == len(record["alt"]) == 1:
        return cigar_observation(read, record)
    if read.is_unmapped or not read.query_sequence or read.reference_name != record["chrom"]:
        return None
    pos, ref, _ = minimal_edit(record)
    start, end = pos - 1, pos - 1 + len(ref)
    aligned = {}
    q, g = 0, read.reference_start
    for op, size in read.cigartuples:
        if op == 3 and g <= end and g + size > start - 1:
            return None
        if op in (0, 7, 8):
            aligned.update((g + i, q + i) for i in range(size))
        if op in (0, 1, 4, 7, 8):
            q += size
        if op in (0, 2, 3, 7, 8):
            g += size
    if start - 1 not in aligned or end not in aligned:
        return None
    first, last = aligned[start - 1] + 1, aligned[end]
    return dict(record=record_digest(read), name=read.query_name,
                allele=read.query_sequence[first:last], query_interval=[first, last])


def direct_evidence(reads, record):
    """Record-level primary evidence; no mate merging or normalization."""
    _, ref, alt = minimal_edit(record)
    groups = Counter()
    observations = []
    for read in reads:
        if read.flag & (4 | 256 | 512 | 1024 | 2048) or read.mapping_quality < 1:
            continue
        observation = direct_observation(read, record)
        if observation is None or "N" in observation["allele"]:
            continue
        group = "ref" if observation["allele"] == ref else "alt" if observation["allele"] == alt else "other"
        groups[group] += 1
        observations.append({k: observation[k] for k in ("record", "name", "allele", "query_interval")})
    return dict(counts={g: groups[g] for g in ("ref", "alt", "other")},
                observations=sorted(observations, key=lambda o: (o["record"], o["query_interval"])))


def canonical_audit(result):
    """Keep behavioral outputs, not incidental transcript contribution order."""
    for protein in result.get("proteins", []) + result.get("uncapped_ranked_proteins", []):
        protein["checks"].sort(key=lambda c: (c.get("transcript_version", ""), c.get("protein_start_1based", 0)))
    return result
