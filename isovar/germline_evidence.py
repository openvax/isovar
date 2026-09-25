"""Reads at matched germline variants on the fragments behind a somatic alt allele.

Isovar collects reads only at its input variants, so on its own it cannot say
that a matched germline variant is on the other copy from a somatic variant.
`germline_read_evidence` finds the germline variants inside the aligned blocks
of each somatic variant's alt reads, never inside a skipped intron, and
collects their reads with the same `ReadCollector`, keeping only reads from
fragments that also carry the somatic alt allele.

Only those fragments say which copy the somatic variant is on: the germline
allele beside it is that copy's allele. A fragment with the somatic reference
allele is not evidence either way, because the germline alt allele also comes
from cells without the somatic variant (normal contamination, other subclones)
and from both copies of a homozygous variant. `IsovarReadPhasing.in_cis`
therefore counts somatic-alt fragments carrying the germline alt (cis) against
those carrying the germline reference (trans).

A germline locus reached only by an unmerged mate is not examined, since the
mate's placement is not among the reads at the somatic locus.
"""

from bisect import bisect_left
from collections import defaultdict
import re

from .logging import get_logger
from .read_evidence import ReadEvidence
from .read_identity import fragment_ids
from .variant_helpers import base0_interval_for_variant

logger = get_logger(__name__)

_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def _reference_blocks(start, cigar):
    """[start, end) reference intervals a placement spans, split at skipped introns."""
    blocks, block_start, position = [], start, start
    for length, operation in _CIGAR.findall(cigar or ""):
        if operation in "MD=X":
            position += int(length)
        elif operation == "N":
            if position > block_start:
                blocks.append((block_start, position))
            position += int(length)
            block_start = position
    if position > block_start:
        blocks.append((block_start, position))
    return blocks


def _covers(block_start, block_end, start, end):
    """Whether a block covers a variant's reference interval (both flanks of an insertion)."""
    if start == end:
        return block_start < start < block_end
    return start < block_end and block_start < end


class _GermlineIndex(object):
    """Germline variants by alignment reference name, sorted by start."""

    def __init__(self, germline_variants, references, read_collector):
        by_contig = defaultdict(list)
        for variant in germline_variants:
            by_contig[variant.contig].append(variant)
        rows_by_reference = defaultdict(list)
        for contig, variants in by_contig.items():
            try:
                reference = read_collector._infer_chromosome_name(contig, references)
            except ValueError as error:
                logger.warning("Germline variants on %s are not examined for phasing: %s", contig, error)
                continue
            if reference is not None:
                rows_by_reference[reference].extend(base0_interval_for_variant(v) + (v,) for v in variants)
        self.by_reference = {}
        for reference, rows in rows_by_reference.items():
            rows.sort(key=lambda row: (row[0], row[1]))
            longest = max(end - start for start, end, _ in rows)
            self.by_reference[reference] = ([row[0] for row in rows], rows, longest)

    def covered(self, reference, blocks):
        """Germline variants inside any of these blocks on one reference."""
        if reference not in self.by_reference:
            return set()
        starts, rows, longest = self.by_reference[reference]
        found = set()
        for block_start, block_end in blocks:
            i = bisect_left(starts, block_start - longest)
            while i < len(rows) and rows[i][0] < block_end:
                start, end, variant = rows[i]
                if _covers(block_start, block_end, start, end):
                    found.add(variant)
                i += 1
        return found


def _on_fragments(evidence, fragments):
    """``evidence`` restricted to reads from these fragments."""
    def keep(reads):
        return [read for read in reads if fragment_ids((read,)) & fragments]
    return ReadEvidence(
        trimmed_base1_start=evidence.trimmed_base1_start, trimmed_ref=evidence.trimmed_ref,
        trimmed_alt=evidence.trimmed_alt, ref_reads=keep(evidence.ref_reads),
        alt_reads=keep(evidence.alt_reads), other_reads=keep(evidence.other_reads))


def germline_read_evidence(results, germline_variants, alignment_file, read_collector):
    """
    Reads at the matched germline variants each result's alt reads cover.

    Parameters
    ----------
    results : list of IsovarResult
    germline_variants : sequence of varcode.Variant
        Variants from a matched normal sample. (`run_isovar` rejects a
        germline variant that is also an input; called directly, such
        variants are skipped, since their own results hold their reads.)
    alignment_file : pysam.AlignmentFile
        The alignments the results were collected from.
    read_collector : ReadCollector
        The collector the results were collected with.

    Returns
    -------
    list of dict
        One per result: each covered germline variant to a ReadEvidence of
        its reads from fragments that also carry the result's alt allele.
        Reads without alignment metadata cover nothing.
    """
    in_run = {result.variant for result in results}
    germline_variants = [v for v in germline_variants if v not in in_run]
    if not germline_variants:
        return [{} for _ in results]
    references = alignment_file.references
    index = _GermlineIndex(germline_variants, set(references), read_collector)
    covered_by_result = []
    for result in results:
        blocks = defaultdict(set)
        for read in result.read_evidence.alt_reads:
            for _, (reference_id, start, cigar, _) in getattr(read, "source_alignments", ()):
                blocks[references[reference_id]].update(_reference_blocks(start, cigar))
        covered = set()
        for reference, reference_blocks in blocks.items():
            covered |= index.covered(reference, sorted(reference_blocks))
        covered_by_result.append(sorted(covered, key=lambda v: (v.contig, v.start, v.ref, v.alt)))
    # Collect each locus once, and release its reads after the last result that needs them.
    remaining = defaultdict(int)
    for covered in covered_by_result:
        for variant in covered:
            remaining[variant] += 1
    collected, evidence = {}, []
    for result, covered in zip(results, covered_by_result):
        fragments = fragment_ids(result.read_evidence.alt_reads)
        kept = {}
        for variant in covered:
            if variant not in collected:
                collected[variant] = read_collector.read_evidence_for_variant(variant, alignment_file)
            kept[variant] = _on_fragments(collected[variant], fragments)
            remaining[variant] -= 1
            if not remaining[variant]:
                del collected[variant]
        evidence.append(kept)
    return evidence
