# Incremental audit result decoding (#496)

Read JSON in fixed-size text chunks, retaining at most the current scalar token
besides the materialized result. Delegate scalar decoding to Python's standard
JSON decoder to preserve numeric, escape, surrogate and duplicate-key behavior.
Decode complete small subtrees already in a buffered chunk with the standard
C decoder; never grow a buffer to finish a container. Disable speculative
subtree decoding when a caller supplies an object hook, preserving exactly
one callback for each object. Keep the fallback grammar iterative and validate
one complete document. Intern only
a bounded set of repeated object keys. Existing gzip/plain file formats and
historical canonical bytes remain unchanged.

Compare expected replay results by reducing parser events directly to the
existing semantic fingerprint. Arrays need a running hash; objects retain only
member keys and child hashes until their keys can be sorted, including last-key
wins for duplicates. Do not materialize the expected reconstruction. Preserve
literal-null expected-result behavior and reject malformed trailing content.

A fully materialized read or resumed API result remains output-sized; this is
not a claim of bounded total memory. Do not change the frozen ONT producer or
its environment. Validate chunk-only and split-token streams, stdlib value and
fingerprint equivalence, malformed/truncated inputs, non-ASCII/surrogate values,
large integers, floating-point rounding, duplicate keys and existing replay.
Measure isolated existing-result decoding and comparison before release.

The optional ijson 3.6.0 probe remains outside the repository environment. Its
C backend replaced an unpaired escaped surrogate with a question mark in the
compatibility probe, so it cannot preserve this reader's existing value contract.
