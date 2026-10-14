# Batch standard-library JSON tokens (#504)

The local feasibility experiment found that joining small encoder tokens before
UTF-8 conversion reduced compressed writing time by about 19% on a complete
18.3 MB ONT reconstruction. Standard-library buffered text output improved it
by about 9%; whole-document `json.dumps` was faster but allocated a serialized
copy proportional to the result size. Five interleaved trials measured writing,
not reconstruction. Synthetic small-field evidence also benefited; large strings
showed little improvement.

Keep `JSONEncoder.iterencode` and the existing canonical options. Accumulate
small tokens up to the existing chunk threshold, then join and encode once.
Flush pending tokens before slicing an oversized token. Reuse this through the
existing hashing and writing callers. Add no serializer, traversal framework,
dependency, public option or cache.

Preserve complete JSON and gzip bytes, checksums and standard-encoder errors.
The pending buffer remains bounded independently of document size. As before,
the standard encoder can materialize one oversized string token; this change
does not add a stricter guarantee for individual huge strings.

Verify surrounding punctuation and Unicode when oversized strings interrupt
small tokens, run the existing audit/checkpoint memory and failure checks, and
compare the full saved ONT output with the released writer. Record benchmark
methodology and measurements in the PR. Run lint, the full suite and four CI
jobs, then deploy the patch release from clean master and verify the live wheel.
