
from collections import defaultdict

nucleotides = ["A", "C", "T", "G"]
nucleotide_to_index = {c: i for (i, c) in enumerate(nucleotides)}
index_to_nucleotide = {i: c for (i, c) in enumerate(nucleotides)}

def group_unique_sequences(
        variant_reads,
        max_prefix_size=None,
        max_suffix_size=None):
    """
    Given a list of VariantRead objects, extracts all unique
    (prefix, suffix) sequence pairs and associate each with a list
    of read names that contained that sequence.
    """
    groups = defaultdict(set)
    for r in variant_reads:
        prefix = r.prefix
        suffix = r.suffix
        if max_prefix_size and len(prefix) > max_prefix_size:
            prefix = prefix[-max_prefix_size:]
        if max_suffix_size and len(suffix) > max_suffix_size:
            suffix = suffix[:max_suffix_size]
        key = (prefix, suffix)
        groups[key].add(r.name)
    return groups

def count_unique_sequences(
        variant_reads,
        max_prefix_size=None,
        max_suffix_size=None):
    """
    Given a list of VariantRead objects, extracts all unique
    (prefix, suffix) sequence pairs and associate each with the number
    of reads that contain that sequence.
    """
    groups = group_unique_sequences(
        variant_reads,
        max_prefix_size=max_prefix_size,
        max_suffix_size=max_suffix_size)
    return {
        seq_pair: len(read_names)
        for (seq_pair, read_names) in groups.items()
    }


def get_variant_nucleotides(variant_reads):
    if len(variant_reads) > 0:
        variant_seq = variant_reads[0].variant
        if not all(r.variant == variant_seq for r in variant_reads):
            raise ValueError(
                ("Cannot call `get_variant_nucleotides` on a collection"
                 " of VariantRead objects spanning multiple variants"))
        return variant_seq
    else:
        raise ValueError("Expected len(variant_reads) > 0")


def make_prefix_suffix_pairs(variant_reads):
    variant_seq = get_variant_nucleotides(variant_reads)
    pairs = [(r.prefix, r.suffix) for r in variant_reads]
    return variant_seq, pairs
