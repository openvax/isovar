
from collections import Counter

nucleotide_to_index = {"A": 0, "C": 1, "G": 2, "T": 3}
index_to_nucleotide = {i: c for (c, i) in nucleotide_to_index.items()}


def unique_counts(tuple_list):
    """
    Given a list of (prefix, suffix) sequence parts,
    return a Counter of distinct sequences
    """
    counter = Counter()
    for x in tuple_list:
        counter[x] += 1
    return counter

def drop_variant_from_partitioned_sequences(partitioned_read_sequences):
    assert len(partitioned_read_sequences) > 0
    variant_seq = partitioned_read_sequences[0][1]
    # get rid of the repeated variant nucleotides and keep only the
    # prefix and suffix
    pairs = [(p, s) for (p, _, s) in partitioned_read_sequences]
    return variant_seq, pairs
