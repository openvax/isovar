"""Loose, positional initiation priors, separate from translation observations."""

from dataclasses import asdict, dataclass
from typing import Optional


# A of ATG is +1; there is no position zero. The outer GCC is optional.
KOZAK_PREFERENCES = {
    -9: "G", -8: "C", -7: "C", -6: "G", -5: "C", -4: "C",
    -3: "AG", -2: "C", -1: "C", 4: "G",
}


@dataclass(frozen=True)
class InitiationContext:
    """Transcript-oriented bases surrounding a proposed initiation codon.

    Missing flanks are empty/short strings, never inferred from protein
    sequence or intronic genomic adjacency. Ambiguous bases remain unknown.
    ``genomic_positions`` are optional 1-based placements in sequence order.
    ``start_offset`` is measured from the supplied transcript sequence; it
    is not a verified cap distance unless ``five_prime_complete`` is True.
    """

    upstream: str
    start_codon: str
    downstream: str
    genomic_positions: tuple = ()
    source: Optional[str] = None
    source_transcript_id: Optional[str] = None
    start_offset: Optional[int] = None
    five_prime_complete: Optional[bool] = None
    unavailable_reason: Optional[str] = None

    def __post_init__(self):
        sequence = self.upstream + self.start_codon + self.downstream
        if (len(self.upstream) > 9 or len(self.downstream) > 6
                or any(b not in "ACGTRYSWKMBDHVN" for b in sequence)
                or (len(self.start_codon) != 3 and not self.unavailable_reason)):
            raise ValueError("Invalid initiation-context sequence")
        if self.unavailable_reason and sequence:
            raise ValueError("Unavailable context must not carry unvalidated sequence")
        if self.genomic_positions and (
                len(self.genomic_positions) != len(sequence)
                or any(type(p) is not int or p < 1 for p in self.genomic_positions)
                or len(set(self.genomic_positions)) != len(self.genomic_positions)):
            raise ValueError("Invalid initiation-context genomic placements")
        if self.start_offset is not None and (
                type(self.start_offset) is not int or self.start_offset < len(self.upstream)):
            raise ValueError("Invalid initiation-context start offset")
        if self.five_prime_complete is not None and type(self.five_prime_complete) is not bool:
            raise ValueError("five_prime_complete must be True, False or None")

    @property
    def relative_positions(self):
        """Kozak numbering for every retained base, without a position zero."""
        return tuple(range(-len(self.upstream), 0)) + tuple(range(1, len(self.start_codon) + len(self.downstream) + 1))

    @property
    def assessment(self):
        """Return observed preferences, missingness and a qualitative prior.

        Match counts are descriptive, not calibrated probabilities or a
        translation threshold. A mismatch never establishes no initiation.
        Non-ATG starts are retained without predicting their efficiency.
        """
        sequence = self.upstream + self.start_codon + self.downstream
        bases = dict(zip(self.relative_positions, sequence))
        positions = []
        for position, preferred in KOZAK_PREFERENCES.items():
            base = bases.get(position)
            positions.append(dict(position=position, base=base, preferred=preferred,
                                  matches=base in preferred if base is not None and base in "ACGT" else None))
        key = [p for p in positions if p["position"] in (-3, 4)]
        observed = sum(p["matches"] is not None for p in key)
        matches = sum(p["matches"] is True for p in key)
        if self.unavailable_reason:
            status, interpretation = "unavailable", "sequence_context_unavailable"
        else:
            status = "complete" if len(self.upstream) == 9 and len(self.downstream) == 6 and all(b in "ACGT" for b in sequence) else "partial"
            interpretation = (
                "both_key_preferences_present" if matches == 2 else
                "one_key_preference_present" if matches == 1 else
                "no_observed_key_preference; initiation_still_possible" if observed else
                "key_positions_unknown; initiation_unresolved")
        codon = self.start_codon
        atg = None if len(codon) != 3 or any(b not in "ACGT" for b in codon) else codon == "ATG"
        boundary = ("at_verified_5prime_end; conventional_scanning_uncertain" if self.five_prime_complete else
                    "at_supplied_sequence_start; true_5prime_end_unverified") if self.start_offset == 0 else None
        return dict(policy="isovar.initiation_context.v1", status=status,
                    display=self.upstream + "[" + codon + "]" + self.downstream if sequence else None,
                    is_atg=atg, kozak_positions=positions,
                    key_positions_observed=observed, key_positions_matching=matches,
                    interpretation=interpretation, sequence_boundary_note=boundary,
                    cap_distance_nt=self.start_offset if self.five_prime_complete is True else None,
                    translation_evidence="sequence_prior_only; not_observed")

    def to_dict(self):
        """Serialize raw sequence/provenance and the reproducible assessment."""
        result = asdict(self)
        result["assessment"] = self.assessment
        return result


def annotate_initiation_context(sequence, start, *, five_prime_complete=None):
    """Extract up to 9 upstream and 6 downstream bases from a supplied RNA.

    Parameters
    ----------
    sequence : str
        Transcript-oriented sequence, accepting RNA U or DNA T and IUPAC
        ambiguity. It need not be a complete mature mRNA.
    start : int
        Zero-based offset of a complete proposed initiation codon. ATG is
        not required; retained non-ATG codons receive no efficiency claim.
    five_prime_complete : bool or None
        Caller evidence that the sequence begins at the biological 5-prime
        end. The default leaves the cap distance unknown.

    Returns
    -------
    InitiationContext
        Observed flanks and a loose positional assessment. Partial contexts
        and ATG at the first three bases are retained, never filtered out.
    """
    sequence = sequence.upper().replace("U", "T")
    if type(start) is not int or not 0 <= start <= len(sequence) - 3:
        raise ValueError("start must identify a complete codon in the sequence")
    return InitiationContext(sequence[max(0, start - 9):start], sequence[start:start + 3],
                             sequence[start + 3:start + 9], start_offset=start,
                             five_prime_complete=five_prime_complete)
