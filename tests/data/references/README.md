# Native mitochondrial fixtures

`mitochondrial.json` is copied unchanged from ORFdata's pinned reference-binding
fixtures (openvax/orfdata, commit 2ac57b6). It includes source URLs, checksums,
native coordinates and protein annotations for human OpenProt 2.2 IP_306399
(Ensembl 106 MT-ND4) and yeast OLI1/Q0130 and Q0092 (Ensembl 115).

OpenProt data are CC-BY-4.0; Ensembl data follow its data-use terms. The source
and evidence metadata remain in the fixture. Tests use synthetic nucleotide
changes on these native contexts; they do not represent new RNA or peptide
observations. MS/translation-event scores are aggregate source counts.
