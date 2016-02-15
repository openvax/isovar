
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
        if max_prefix_size:
            prefix = prefix[-max_prefix_size:]
        if max_suffix_size:
            suffix = suffix[-max_suffix_size:]
        key = (r.prefix, r.suffix)
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

def make_prefix_suffix_pairs(variant_reads):
    assert len(variant_reads) > 0
    variant_seq = variant_reads[0].variant
    pairs = [(r.prefix, r.suffix) for r in variant_reads]
    return variant_seq, pairs
