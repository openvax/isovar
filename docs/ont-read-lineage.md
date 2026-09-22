# ONT read lineage

Implementation spec for [#331](https://github.com/openvax/isovar/issues/331):

- Keep segment, fragment, alignment compatibility and reconstruction thresholds
  unchanged. Add lineage evidence to SV junction support and to already verified
  full-interval ORF witnesses. Partial children cannot supply a complete witness.
- Recognize Dorado only with an ONT read group and a documented program chain.
  Follow explicit record/RG program pointers; otherwise require every header
  program chain to originate from Dorado basecalling (`PN:dorado` and a `CL`
  running `basecaller` or `duplex`). Dorado alignment alone is insufficient.
  Never infer a producer from tag names.
- Group simplex split children by their explicit `pi` parent, within a read
  group and the current input. Follow visible parent relationships and reject
  cycles, conflicting tags and unresolved parents. Retain `dx`/`pi`/`sp` evidence.
- Report duplex parents/consensuses as unresolved: `dx` gives an output class,
  not the parent-to-consensus mapping. Do not parse QNAMEs to invent that mapping.
- Cache tag interpretation per segment; do not walk sequences or decode arrays.
  Test producer ambiguity, namespaces, split families, duplex mixtures,
  alternative alignments and complete ORF witnesses. Existing Sid records remain
  unchanged; they do not carry the Dorado lineage tags needed for these cases.

The reported groups describe **sequencing signal ancestry**, not independent RNA
molecules. Dorado can split concatenated molecules in one input signal, and PCR
copies can produce separate signals. Read groups remain separate even when their
`LB` strings match; merging lanes or technical reprocessing needs the explicit
library/sample contract in [#226](https://github.com/openvax/isovar/issues/226).

Primary definitions: [Dorado SAM tags](https://software-docs.nanoporetech.com/dorado/latest/basecaller/sam_spec/),
[read splitting](https://software-docs.nanoporetech.com/dorado/latest/basecaller/read_splitting/),
[duplex output](https://software-docs.nanoporetech.com/dorado/latest/basecaller/duplex/),
and [SAM program/read-group fields](https://samtools.github.io/hts-specs/SAMv1.pdf).
