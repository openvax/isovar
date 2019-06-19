from isovar.variant_sequence import VariantSequence
from isovar.assembly import iterative_overlap_assembly

import random
import time

n_reads = 2000
read_length = 150
full_sequence_length = 190
variant_pos = 100

full_sequence = "".join(random.choices("ACTG", k=full_sequence_length))

sequences = []
for i in range(n_reads):
    start_pos = random.randint(0, full_sequence_length - read_length)
    prefix = full_sequence[start_pos:variant_pos]
    alt = full_sequence[variant_pos]
    suffix = full_sequence[variant_pos:start_pos + read_length]
    seq = VariantSequence(
        prefix=prefix,
        alt=alt,
        suffix=suffix,
        reads={"read_%d" % i})
    sequences.append(seq)


t0 = time.time()
results = iterative_overlap_assembly(sequences)
t1 = time.time()
print("Assembled %d sequences, time = %0.2f" % (
    len(results),
    t1 - t0))
