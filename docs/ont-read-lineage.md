# ONT read lineage

Oxford Nanopore's Dorado basecaller can turn one raw signal into several reads.
It splits concatenated molecules into child reads, and duplex calling adds a
consensus read made from two strands. These reads have different names but
share their signal, so counting them as separate reads can overstate support.

For SV junctions, and for reads covering a whole ORF, Isovar reports this signal
lineage alongside the read and fragment counts, never in place of them. Lineage
does not change which alignments are compatible or any reconstruction threshold.
A partial child read cannot stand in for a read covering the whole ORF.

## When lineage is used

- **Producer.** Lineage tags are read only for an ONT read group whose header
  shows a Dorado program chain. Isovar follows explicit record or read-group
  `PG` pointers. Without them, every program chain in the header must start with
  Dorado basecalling: `PN:dorado` with a `CL` running `basecaller` or `duplex`.
  Alignment by Dorado alone is not enough, and the producer is never guessed from
  tag names.
- **Split reads.** Simplex child reads are grouped by their explicit `pi` parent,
  within one read group of one input. Cycles, conflicting tags and unresolved
  parents are reported as unresolved, and the `dx`/`pi`/`sp` tags are kept.
- **Duplex reads.** Duplex parents and consensus reads are reported as
  unresolved. `dx` gives the read's class but not which parents made which
  consensus, and read names are not parsed to guess it.

## What lineage does not tell you

The groups describe **sequencing signal ancestry**, not independent RNA
molecules. One signal can contain several molecules that Dorado splits, and PCR
copies of one molecule produce separate signals. Read groups stay separate even
when their `LB` strings match. Merging lanes or reprocessed data needs the
explicit library and sample contract in
[#226](https://github.com/openvax/isovar/issues/226).

The packaged Sid test reads do not carry Dorado lineage tags. Synthetic tests
cover producer ambiguity, split families, duplex mixtures and alternative
alignments.

Primary definitions: [Dorado SAM tags](https://software-docs.nanoporetech.com/dorado/latest/basecaller/sam_spec/),
[read splitting](https://software-docs.nanoporetech.com/dorado/latest/basecaller/read_splitting/),
[duplex output](https://software-docs.nanoporetech.com/dorado/latest/basecaller/duplex/),
and [SAM program and read group fields](https://samtools.github.io/hts-specs/SAMv1.pdf).
