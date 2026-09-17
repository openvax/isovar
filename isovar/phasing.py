# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from collections import defaultdict, Counter
from itertools import combinations

from .default_parameters import MIN_SHARED_FRAGMENTS_FOR_PHASING
from .phase_group import PhaseGroup
from .read_identity import alignment_constraints, fragment_ids
from .chimeric_alignment import compatible_phasing_alignments
from .transcript_edit_helpers import transcript_assembly_edit_sort_key


def _variant_sort_key(variant):
    return (
        variant.contig,
        variant.start,
        variant.ref,
        variant.alt,
    )


def create_read_names_to_variants_dict(variant_to_read_names_dict):
    """
    Invert a variant -> read-name mapping into read-name -> variants.
    """
    read_names_to_variants = defaultdict(set)

    for variant, read_names in variant_to_read_names_dict.items():
        for read_name in read_names:
            read_names_to_variants[read_name].add(variant)
    return read_names_to_variants


def create_variant_to_alt_read_names_dict(isovar_results):
    """
    Create dictionary from variant to names of alt reads supporting
    that variant in an IsovarResult.

    Parameters
    ----------
    isovar_results : list of IsovarResult

    Returns
    -------
    Dictionary from varcode.Variant to set(str) of read names

    """
    return {
        isovar_result.variant: set(isovar_result.alt_read_names)
        for isovar_result in isovar_results
    }


def create_variant_to_protein_sequence_read_names_dict(isovar_results):
    """
    Create dictionary from variant to names of alt reads used to create
    mutant protein sequence in an IsovarResult.

    Parameters
    ----------
    isovar_results : list of IsovarResult

    Returns
    -------
    Dictionary from varcode.Variant to set(str) of read names
    """
    variant_to_read_names = {}
    for isovar_result in isovar_results:
        if isovar_result.has_mutant_protein_sequence_from_rna:
            protein_sequence = isovar_result.top_protein_sequence
            read_names = set(
                protein_sequence.read_names_supporting_protein_sequence)
        else:
            read_names = set()
        variant_to_read_names[isovar_result.variant] = read_names
    return variant_to_read_names


def create_variant_to_top_protein_sequence_dict(isovar_results):
    """
    Create dictionary from variant to its top translated protein sequence.

    Variants without an assembled protein sequence are mapped to None.
    """
    return {
        isovar_result.variant: (
            isovar_result.top_protein_sequence
            if isovar_result.has_mutant_protein_sequence_from_rna
            else None
        )
        for isovar_result in isovar_results
    }


def compute_phasing_counts(variant_to_read_names_dict):
    """

    Parameters
    ----------
    variants_to_read_names : dict
        Dictionary mapping varcode.Variant to a set of hashable fragment IDs.
        Plain read names remain supported for caller-provided legacy mappings.

    Returns
    -------
    Dictionary from variant to Counter(Variant)
    """
    support, _ = _phasing_support(variant_to_read_names_dict)
    return defaultdict(Counter, {
        variant: Counter({other: len(ids) for other, ids in neighbors.items()})
        for variant, neighbors in support.items()
    })

def threshold_phased_variant_counts(counts_dict, min_count):
    """
    Choose set of phased variants by keeping any variants with associated
    counts greater than or equal the given threshold.

    Parameters
    ----------
    counts_dict : variant -> int dict

    min_count : int

    Returns
    -------
    set of varcode.Variant
    """
    return {
        variant
        for (variant, count)
        in counts_dict.items()
        if count >= min_count
    }


def create_phase_groups(
        variant_to_read_names_dict,
        min_shared_fragments_for_phasing,
        variant_to_top_protein_sequence_dict=None,
        read_names_by_id=None):
    """
    Group variants into connected components of the phasing graph.

    If top translated protein sequences are provided then each resulting
    PhaseGroup is also annotated with directly observed cDNA, protein, and
    transcript metadata from those assemblies.

    ``read_names_by_id`` optionally maps scoped fragment IDs to display names.
    Names are converted only after constructing the graph; public PhaseGroup
    fields retain their existing string-based representation.

    Returns
    -------
    dict
        Mapping from variant to PhaseGroup. Variants without phased partners are
        omitted.
    """
    return _phase_annotations(
        variant_to_read_names_dict,
        min_shared_fragments_for_phasing,
        variant_to_top_protein_sequence_dict,
        read_names_by_id)[1]


def _phasing_support(variant_to_reads):
    """Build pairwise edges from compatible observations, once per fragment.

    Alternative observations are hypotheses, not simultaneous constraints.
    A pair needs at least one jointly compatible choice. Reuse the assembly
    rule for linear placements and complementary mates. Different placements
    of one segment require reciprocal, compatible supplementary-path evidence.

    Metadata-free reads/IDs retain legacy name-only semantics. Public names
    are kept separately and never substituted for scoped evidence IDs.
    """
    observations = defaultdict(lambda: defaultdict(list))
    names_by_id = {}
    for variant, reads in variant_to_reads.items():
        for read in reads:
            constraints = alignment_constraints((read,))
            for fragment in fragment_ids((read,)):
                observations[fragment][variant].append((constraints, read))
                names_by_id[fragment] = getattr(read, "name", read)
    support = defaultdict(dict)
    for fragment, by_variant in observations.items():
        placements = defaultdict(set)
        for reads in by_variant.values():
            for constraints, _ in reads:
                for segment, placement in constraints.items():
                    placements[segment].add(placement)
        for (variant, reads), (other, other_reads) in combinations(by_variant.items(), 2):
            if any(compatible_phasing_alignments(constraints, read, other_read, placements)
                   for constraints, read in reads for _, other_read in other_reads):
                shared = support[variant].setdefault(other, set())
                shared.add(fragment)
                support[other][variant] = shared
    return support, names_by_id


def _phase_annotations(
        variant_to_reads,
        min_shared_fragments_for_phasing,
        variant_to_top_protein_sequence_dict=None,
        read_names_by_id=None):
    """Derive neighbors and groups from the same validated fragment edges."""
    support, names_by_id = _phasing_support(variant_to_reads)
    if read_names_by_id is not None:
        names_by_id.update(read_names_by_id)
    phased_neighbors = {
        variant: {other for other, ids in support[variant].items()
                  if len(ids) >= min_shared_fragments_for_phasing}
        for variant in variant_to_reads
    }

    visited = set()
    variant_to_phase_group = {}
    for variant in sorted(variant_to_reads, key=_variant_sort_key):
        if variant in visited:
            continue

        component = set()
        pending = [variant]
        while pending:
            current_variant = pending.pop()
            if current_variant in component:
                continue
            component.add(current_variant)
            visited.add(current_variant)
            pending.extend(phased_neighbors.get(current_variant, set()) - component)

        if len(component) <= 1:
            continue

        supporting_read_names = {
            names_by_id[fragment]
            for grouped_variant in component
            for other in phased_neighbors[grouped_variant]
            for fragment in support[grouped_variant][other]
        }

        cdna_sequences = set()
        mutant_protein_sequences = set()
        transcript_ids = set()
        transcript_names = set()
        known_somatic_transcript_edits = set()
        known_germline_transcript_edits = set()
        unexplained_transcript_edits = set()
        for grouped_variant in component:
            protein_sequence = (variant_to_top_protein_sequence_dict or {}).get(grouped_variant)
            if protein_sequence is None:
                continue
            if hasattr(protein_sequence, "_transcript_assembly_edits_by_category"):
                categorized_edits = protein_sequence._transcript_assembly_edits_by_category()
            else:
                categorized_edits = {
                    "known_somatic": getattr(protein_sequence, "known_somatic_transcript_edits", ()),
                    "known_germline": getattr(protein_sequence, "known_germline_transcript_edits", ()),
                    "unexplained": getattr(protein_sequence, "unexplained_transcript_edits", ()),
                }
            cdna_sequences.update(protein_sequence.cdna_sequences)
            mutant_protein_sequences.add(protein_sequence.amino_acids)
            transcript_ids.update(protein_sequence.transcript_ids)
            transcript_names.update(protein_sequence.transcript_names)
            known_somatic_transcript_edits.update(categorized_edits["known_somatic"])
            known_germline_transcript_edits.update(categorized_edits["known_germline"])
            unexplained_transcript_edits.update(categorized_edits["unexplained"])

        phase_group = PhaseGroup(
            somatic_variants=tuple(sorted(component, key=_variant_sort_key)),
            germline_variants=(),
            supporting_read_names=supporting_read_names,
            cdna_sequences=tuple(sorted(cdna_sequences)),
            mutant_protein_sequences=tuple(sorted(mutant_protein_sequences)),
            transcript_ids=tuple(sorted(transcript_ids)),
            transcript_names=tuple(sorted(transcript_names)),
            known_somatic_transcript_edits=tuple(sorted(
                known_somatic_transcript_edits,
                key=transcript_assembly_edit_sort_key,
            )),
            known_germline_transcript_edits=tuple(sorted(
                known_germline_transcript_edits,
                key=transcript_assembly_edit_sort_key,
            )),
            unexplained_transcript_edits=tuple(sorted(
                unexplained_transcript_edits,
                key=transcript_assembly_edit_sort_key,
            )),
        )
        for grouped_variant in component:
            variant_to_phase_group[grouped_variant] = phase_group
    return phased_neighbors, variant_to_phase_group


def _variant_reads(isovar_results, protein=False):
    """Retain placement evidence; accept names-only legacy result objects."""
    by_variant = {}
    for result in isovar_results:
        if protein:
            if result.has_mutant_protein_sequence_from_rna:
                sequence = result.top_protein_sequence
                reads = getattr(sequence, "supporting_reads", None)
                if reads is None:
                    reads = sequence.read_names_supporting_protein_sequence
            else:
                reads = ()
        else:
            reads = getattr(result, "alt_reads", None)
            if reads is None:
                reads = result.alt_read_names
        by_variant[result.variant] = reads
    return by_variant


def annotate_phased_variants(
        unphased_isovar_results,
        min_shared_fragments_for_phasing=MIN_SHARED_FRAGMENTS_FOR_PHASING):
    """
    Annotate IsovarResult objects with phasing information. Phasing
    is determined by looking at RNA fragments used for assembled protein
    sequences and counting the number of shared fragments with compatible
    placements, scoped by SAM read group. Public read-name sets remain display
    names, not evidence IDs. Legacy caller-created objects without alignment
    metadata retain name-only phasing;
    they are not equated with scoped collected fragments.

    Parameters
    ----------
    unphased_isovar_results : list of IsovarResult

    min_shared_fragments_for_phasing : int

    Returns
    -------
    list of IsovarResult
    """

    updates = {result.variant: {} for result in unphased_isovar_results}
    for source in ("supporting_reads", "protein_sequence"):
        protein = source == "protein_sequence"
        neighbors, groups = _phase_annotations(
            _variant_reads(unphased_isovar_results, protein=protein),
            min_shared_fragments_for_phasing,
            create_variant_to_top_protein_sequence_dict(unphased_isovar_results) if protein else None)
        for variant, fields in updates.items():
            fields["phased_variants_in_" + source] = neighbors[variant]
            fields["phase_group_from_" + source] = groups.get(variant)
    return [result.clone_with_updates(**updates[result.variant])
            for result in unphased_isovar_results]
