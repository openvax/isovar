# Recovery event report memory

Related: #495 and #470.

The event report currently retains reconstructed paths from every view, including large support and exploratory ORF payloads that it never uses. Retain only sequence, blocks, end reasons and unplaced intervals for report assembly. Release each full reconstruction before loading another. Keep input checksum validation, path sequence/ORF agreement, coordinates and published values unchanged.

Verify with a two-view lifetime regression and a mismatched ORF sequence rejection. This changes report generation only; the frozen reconstruction producer and its scientific settings remain unchanged.
