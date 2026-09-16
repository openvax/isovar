# Figure notes

Reference: GRCh38-osteosarc-vaccine-cohort-ensembl87-7be1ab94b01efb6b.

Orange marks the mutation or deletion boundary; blue denotes assembly on, gray assembly off. Protein peptide counts are mutation-overlapping windows of the displayed length.

Reconstructed cDNA is one actual RNA-derived sequence producing the displayed protein. Other reconstructions can produce the same protein; their reads are not combined into this track. The read-overlap panel shows the first mode with a protein, normally assembly on.

Read objects are post-mate-merge observations, not independent molecules. Spanning read objects cover the entire displayed reconstructed cDNA; spans are clipped to that sequence. Identical spans are grouped, with multiplicity marked only when greater than one. Coverage uses every supporting object, including groups hidden by the display limit. Protein template counts include all translations contributing that protein; cDNA counts refer only to the displayed reconstruction. Junction counts are retained CIGAR N observations.

The title uses normalized 1-based variant coordinates. Genomic tracks use forward-strand 0-based, half-open coordinates, with compressed genomic gaps marked //; cDNA and protein offsets are transcript-oriented. Angled gray connectors join adjacent annotated exon boundaries; black connectors on the RNA junction row are observed splice junctions. Transcript names come from the same annotation as the ENST IDs. Models do not establish a unique isoform.

These are candidates before result-level filters, not independent validation or a clinical recommendation. Complete sequences, settings, support and limitations are in evidence.json.

## Varcode baseline

These tracks predict the reference transcript plus the nominated variant only, not observed expression or a patient haplotype. The display is anchored to the genomic edit's codon; Varcode effect labels may use a different repeat-normalized protein boundary. Equal local predicted sequences are grouped, not ranked. No predicted residues fill missing Isovar context. Reference tracks do not mark a terminal stop; their ends may simply be the display boundary. Black transcript boxes contribute to a displayed RNA protein; gray boxes are prediction-only candidates. Compatibility can use read alignment outside the displayed retained cDNA; this local track is not a full explanation of every model exclusion.

| Track | Transcript (name) | Varcode effect |
| --- | --- | --- |
| Varcode 1 | ENST00000349792 (PIP5K1A-001) | p.G461fs |
| Varcode 1 | ENST00000368888 (PIP5K1A-003) | p.G474fs |
| Varcode 1 | ENST00000409426 (PIP5K1A-015) | p.G462fs |
| Varcode 1 | ENST00000441902 (PIP5K1A-016) | p.G434fs |
