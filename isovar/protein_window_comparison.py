"""Compare local protein windows using their declared focal nucleotide anchor.

An in-frame focal gap follows the allele length difference; a frameshift tail
stays boundary-aligned. This is not a homology alignment. Unavailable reference
residues and ambiguous amino acids cannot establish sequence identity.
"""

from collections.abc import Sequence
from typing import Literal, TypedDict


ComparisonStatus = Literal["reference_identical", "amino_acid_divergent", "unresolved"]


class ProteinColumn(TypedDict):
    reference: str
    protein: str
    position: int | None
    reference_position: int | None
    focal_deletion: bool


class ProteinComparison(TypedDict):
    status: ComparisonStatus
    changed_residue_count: int | None
    substituted_positions: list[int]
    inserted_positions: list[int]
    deleted_reference_positions: list[int]
    unresolved_alignment_columns: list[int]
    stop_changed: bool | None


def reference_window_stop(protein: str, ends_with_stop: bool, limit: int,
                          observed_stop: bool) -> bool | None:
    """Locate a comparator's stop within a capped window, preserving missing bases.

    When the observed window ends in a stop, a comparator ending at that same
    residue without its next codon cannot establish either stop retention or
    stop gain. A longer translated comparator establishes that codon is sense.
    """
    if ends_with_stop and len(protein) <= limit:
        return True
    if len(protein) > limit or not observed_stop:
        return False
    return None


def aligned_protein_columns(
        reference: str, protein: str, interval: Sequence[int], *,
        allele_length_delta: int = 0, frameshift: bool = False,
        stop: bool = False, reference_stop: bool = False) -> list[ProteinColumn]:
    """Position residues and focal gaps without realigning repeated sequences.

    ``interval`` is the observed protein's zero-based half-open focal interval;
    ``allele_length_delta`` is reference minus observed nucleotide length.
    Missing reference sequence is rendered as a gap with an out-of-range
    reference position, distinct from a justified insertion (position None).
    """
    begin, end = interval
    if not 0 <= begin <= end <= len(protein):
        raise ValueError("Invalid protein mutation interval")
    delta = allele_length_delta // 3 if not frameshift and allele_length_delta % 3 == 0 else 0
    reference_length = len(reference)
    columns: list[ProteinColumn] = []

    def add(reference_position: int | None, position: int | None, deletion: bool = False) -> None:
        columns.append(dict(
            reference=reference[reference_position]
                if reference_position is not None and 0 <= reference_position < reference_length else "-",
            protein=protein[position] if position is not None and position < len(protein) else "-",
            position=position, reference_position=reference_position, focal_deletion=deletion))

    for position in range(begin):
        add(position, position)
    reference_end = max(begin, end + delta)
    for offset in range(max(reference_end - begin, end - begin)):
        ref_position = begin + offset if begin + offset < reference_end else None
        observed_position = begin + offset if begin + offset < end else None
        add(ref_position, observed_position, deletion=observed_position is None and ref_position is not None)
    for position in range(end, len(protein)):
        add(position + delta, position)
    if stop:
        ref_position = len(protein) + delta
        columns.append(dict(
            reference="*" if reference_stop and ref_position == reference_length else
                reference[ref_position] if 0 <= ref_position < reference_length else "-",
            protein="*", position=len(protein), reference_position=ref_position, focal_deletion=False))
    return columns


def compare_protein_windows(
        reference: str | None, protein: str, interval: Sequence[int], *,
        allele_length_delta: int = 0, frameshift: bool = False,
        stop: bool = False, reference_stop: bool | None = None) -> ProteinComparison:
    """Classify a reference comparison independently of a marked mutation span.

    Changed-residue count includes substitutions, inserted residues and deleted
    reference residues. Stop changes are separate. The count is unknown if any
    compared residue is unavailable, even when another verified change already
    establishes divergence. ``reference_stop`` describes the compared window,
    not a stop outside it. RNA sequence comparisons do not demonstrate translation.
    """
    substitutions: list[int] = []
    insertions: list[int] = []
    deletions: list[int] = []
    unresolved: list[int] = []
    if reference is None:
        return dict(status="unresolved", changed_residue_count=None, substituted_positions=[],
                    inserted_positions=[], deleted_reference_positions=[],
                    unresolved_alignment_columns=list(range(len(protein))), stop_changed=None)
    columns = aligned_protein_columns(reference, protein, interval,
                                      allele_length_delta=allele_length_delta, frameshift=frameshift)
    ambiguous = set("XBZJ?-")
    for column_index, column in enumerate(columns):
        observed, expected = column["protein"], column["reference"]
        position, ref_position = column["position"], column["reference_position"]
        if position is None:
            if expected in ambiguous:
                unresolved.append(column_index)
            elif ref_position is not None:
                deletions.append(ref_position)
        elif observed in ambiguous or (ref_position is not None and expected in ambiguous):
            unresolved.append(column_index)
        elif ref_position is None:
            insertions.append(position)
        elif expected != observed:
            substitutions.append(position)
    delta = allele_length_delta // 3 if not frameshift and allele_length_delta % 3 == 0 else 0
    stop_changed = (None if reference_stop is None else
                    stop != reference_stop or bool(stop and reference_stop and
                                                    len(reference) != len(protein) + delta))
    changed = len(substitutions) + len(insertions) + len(deletions)
    status: ComparisonStatus = ("amino_acid_divergent" if changed or stop_changed else
                                "unresolved" if unresolved or stop_changed is None else "reference_identical")
    return dict(status=status, changed_residue_count=None if unresolved else changed,
                substituted_positions=substitutions, inserted_positions=insertions,
                deleted_reference_positions=deletions, unresolved_alignment_columns=unresolved,
                stop_changed=stop_changed)
