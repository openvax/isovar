# ONT read lineage

SV junction support and verified full-interval ORF witnesses report Oxford
Nanopore signal lineage alongside, never instead of, segment and fragment counts.
Alignment compatibility and reconstruction thresholds are unaffected, and a
partial child cannot supply a complete witness.

- **Producer.** Lineage tags are interpreted only for an ONT read group with a
  documented Dorado program chain. Isovar follows explicit record or read-group
  `PG` pointers; otherwise every header program chain must originate from Dorado
  basecalling (`PN:dorado` with a `CL` running `basecaller` or `duplex`). Dorado
  alignment alone is insufficient, and a producer is never inferred from tag names.
- **Split reads.** Simplex split children are grouped by their explicit `pi`
  parent, within one read group and the current input. Cycles, conflicting tags
  and unresolved parents are reported as unresolved, and the `dx`/`pi`/`sp`
  evidence is retained.
- **Duplex.** Duplex parents and consensuses are reported as unresolved: `dx`
  gives an output class, not the parent-to-consensus mapping, and QNAMEs are not
  parsed to invent one.

The packaged Sid records do not carry Dorado lineage tags; synthetic tests cover
producer ambiguity, split families, duplex mixtures and alternative alignments.

The reported groups describe **sequencing signal ancestry**, not independent RNA
molecules. Dorado can split concatenated molecules in one input signal, and PCR
copies can produce separate signals. Read groups remain separate even when their
`LB` strings match; merging lanes or technical reprocessing needs the explicit
library/sample contract in [#226](https://github.com/openvax/isovar/issues/226).

Primary definitions: [Dorado SAM tags](https://software-docs.nanoporetech.com/dorado/latest/basecaller/sam_spec/),
[read splitting](https://software-docs.nanoporetech.com/dorado/latest/basecaller/read_splitting/),
[duplex output](https://software-docs.nanoporetech.com/dorado/latest/basecaller/duplex/),
and [SAM program/read-group fields](https://samtools.github.io/hts-specs/SAMv1.pdf).
