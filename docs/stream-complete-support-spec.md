# Streaming complete-input support (#482)

Freeze discovery and index candidate junctions before scanning original templates.
Build observations in batches bounded by observed bases, consume each batch into
junction counts, linked intervals and exact interval witnesses, then release its
per-base geometry. Keep the existing bounded geometry cache; never keep an
input-wide coordinate interning table. One indivisible template may exceed the
batch budget. No public tuning option or output schema change is needed.

Retain original identities, SAM payloads, qualities, query intervals, alternative
placements, missingness, visible mates and producer-aware ancestry. Resolve
counts from original identities after streaming, so shared cells and parent
signals cannot be double-counted across batches. Direct-junction observations
remain distinct from sequence-exact full-interval witnesses. Discovery budgets,
biological interpretation and all deterministic output ordering stay unchanged.

Validate whole-result equivalence against 1.52.0 on synthetic edge cases and the
four pinned public benchmarks. Exercise tiny batch limits to force boundaries.
Measure transient geometry separately from retained output on synthetic depth
scaling, with explicit runtime, memory and scratch budgets. Run required lint,
full tests, installed-package checks and CI before merging and releasing.

After this release, address #470 in a separate PR/attempt: verify preserved input
and reference identities, declare runtime/scratch budgets before both ONT IGKC
orientations, preserve earlier failures, independently validate original-input
support, and update coverage, event reporting and portable replay provenance.
Unknown outcomes remain unknown until those checks finish.
