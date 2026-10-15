# CTA reconstruction and abundance validation

The initial #465 implementation reuses the coordinate-aware RNA observation,
exact-overlap assembly and disk-backed full-input indexing from SV RNA. It adds
no general transcript quantifier: Salmon estimates abundance against the caller's
full competing transcriptome after target-reference replacement.

## Known mixtures

`examples/benchmark_antigen_isoforms.py` creates two synthetic transcripts with
an alternative cassette exon, shared sequence, and lengths of 450 and 360 nt.
For each of two samples it generates 2,500 paired fragments (100-nt mates,
approximately 200-nt inserts). Sampling probabilities account for effective
length, so fragment proportions are not confused with transcript proportions.
The BAM and FASTQs describe the same original fragments.

| Sample | True T1 transcript fraction | Salmon estimate | Dominant transcript | Accounted fragments |
|---|---:|---:|---|---:|
| A | 0.800 | 0.7982 | T1 | 2,500 / 2,500 |
| B | 0.200 | 0.2056 | T2 | 2,500 / 2,500 |

The full transcript sequences match the independently specified inputs. No
fragment is unassigned or assigned only to a partial sequence in these mixtures.
The integration uses the current stable Salmon 2.8.0 release, downloaded into
the worktree and checked against its GitHub release-asset SHA256. Neither the
system Salmon installation nor the shared Python environment was changed.

The machine-readable record in
[`data/cta-isoforms-experiment.json`](data/cta-isoforms-experiment.json) pins
the script, source modules, Salmon executable archive and input FASTQs/FASTA.
Reproduce from the repository root with an installed Salmon executable:

```sh
PYTHONPATH=. python examples/benchmark_antigen_isoforms.py \
    --salmon /path/to/salmon --output-dir /new/experiment-directory
```

This is a small functional experiment, not a patient result, a real CTA
expression survey, or a general performance benchmark. Real multi-sample CTA
and viral validation remains in #465. No tumor/cell identity or antigen
presentation is inferred from these sequences.

## Regression checks

Offline tests cover both RNA orientations, competing splice isoforms, an
observed new splice, sequence alleles, partial sequence, complete CDS without
UTRs, start loss, missing coding annotation, identical proteins from different
transcripts, related loci, alternative placements, mate/read-group identity,
duplicate payloads, missing qualities and retained barcode labels. Fragment
classes conserve input support; a discovery limit still counts later reads.

An arbitrary 99:1 abundance allocation over indistinguishable RNA models
remains unresolved. Missing competing placements, incomplete discovery and
unresolved partial alleles cannot establish a dominant isoform. Separate
declared samples cannot be pooled by the reconstruction API. Protein exports
omit uncovered reference CDS; protein abundance reports expose the reference
hypothesis contribution separately.

PyEnsembl tests cover genomic exon extraction, CDS offsets, GENCODE-style FASTA
headers, separate Ensembl transcript-version attributes, the normal transcript
sequence API, and BAM/FASTA checksum mismatches despite equal contig lengths.
Reference caches are selected by annotation content and reused across samples.
A regression changes the GTF at the same path and verifies that the new models
are indexed; the underlying stale-index behavior is tracked in
[PyEnsembl #462](https://github.com/openvax/pyensembl/issues/462).
The complete GDC dictionary was extracted from its official GATK-index archive,
verified with published MD5 `f64be73587a7f376c0d8353f1636dca7`. It has 2,779
contigs; the packaged metadata retains archive and dictionary SHA256s.

The GDC genome itself is not bundled or automatically downloaded. The reference
adapter is exercised on small offline DNA/GTF fixtures and the real published
dictionary; this experiment did not align a cohort against the full GDC FASTA.

Repository release dependencies were checked against live PyPI. ORFdata was
updated to 0.27.1 and Varcode to 10.15.1 in that environment; `pip check` passed.
Historical audit and reference pins were not changed.
