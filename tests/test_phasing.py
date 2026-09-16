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

import pytest
from varcode import Variant
from varcode.mutant_transcript import TranscriptEdit

from isovar.allele_read import AlleleRead
from isovar.isovar_result import IsovarResult
from isovar.phasing import (
    annotate_phased_variants,
    create_phase_groups,
    create_variant_to_alt_read_names_dict,
    create_variant_to_protein_sequence_read_names_dict,
)
from isovar.read_evidence import ReadEvidence
from isovar.transcript_assembly_edit import TranscriptAssemblyEdit

from .common import eq_


class DummyProteinSequence(object):
    def __init__(
            self,
            read_names,
            amino_acids,
            cdna_sequences,
            transcript_ids,
            transcript_names,
            known_somatic_transcript_edits=(),
            known_germline_transcript_edits=(),
            unexplained_transcript_edits=()):
        self.read_names_supporting_protein_sequence = set(read_names)
        self.amino_acids = amino_acids
        self.cdna_sequences = set(cdna_sequences)
        self.transcript_ids = list(transcript_ids)
        self.transcript_names = list(transcript_names)
        self.known_somatic_transcript_edits = tuple(known_somatic_transcript_edits)
        self.known_germline_transcript_edits = tuple(
            known_germline_transcript_edits)
        self.unexplained_transcript_edits = tuple(unexplained_transcript_edits)


def make_transcript_assembly_edit(
        transcript_id,
        transcript_name,
        cdna_start,
        cdna_end,
        alt_bases,
        source_variant=None):
    return TranscriptAssemblyEdit(
        transcript_id=transcript_id,
        transcript_name=transcript_name,
        edit=TranscriptEdit(
            cdna_start=cdna_start,
            cdna_end=cdna_end,
            alt_bases=alt_bases,
            source_variant=source_variant,
        ),
    )


def make_isovar_result(
        variant,
        alt_read_names,
        protein_read_names,
        amino_acids=None,
        cdna_sequences=(),
        transcript_ids=(),
        transcript_names=(),
        known_somatic_transcript_edits=(),
        known_germline_transcript_edits=(),
        unexplained_transcript_edits=()):
    read_evidence = ReadEvidence(
        trimmed_base1_start=variant.start,
        trimmed_ref=variant.ref,
        trimmed_alt=variant.alt,
        ref_reads=[],
        alt_reads=[
            AlleleRead(prefix="A", allele=variant.alt, suffix="T", name=read_name)
            for read_name in alt_read_names
        ],
        other_reads=[],
    )
    return IsovarResult(
        variant=variant,
        read_evidence=read_evidence,
        predicted_effect=None,
        sorted_protein_sequences=[DummyProteinSequence(
            read_names=protein_read_names,
            amino_acids=amino_acids if amino_acids is not None else variant.alt,
            cdna_sequences=cdna_sequences,
            transcript_ids=transcript_ids,
            transcript_names=transcript_names,
            known_somatic_transcript_edits=known_somatic_transcript_edits,
            known_germline_transcript_edits=known_germline_transcript_edits,
            unexplained_transcript_edits=unexplained_transcript_edits,
        )],
    )


def test_annotate_phased_variants_creates_explicit_phase_groups():
    v1 = Variant("1", 10, "A", "C", normalize_contig_names=False)
    v2 = Variant("1", 11, "G", "T", normalize_contig_names=False)
    v3 = Variant("1", 12, "T", "G", normalize_contig_names=False)
    germline = Variant("1", 9, "A", "G", normalize_contig_names=False)

    somatic_edit_v1 = make_transcript_assembly_edit("tx-1", "TX1", 10, 11, "C", v1)
    somatic_edit_v2 = make_transcript_assembly_edit("tx-2", "TX2", 20, 21, "T", v2)
    somatic_edit_v3 = make_transcript_assembly_edit("tx-2", "TX2", 30, 31, "G", v3)
    germline_edit = make_transcript_assembly_edit(
        "tx-shared", "TX_SHARED", 7, 8, "G", germline)
    unexplained_edit_v1 = make_transcript_assembly_edit("tx-1", "TX1", 4, 5, "G")
    unexplained_edit_v2 = make_transcript_assembly_edit(
        "tx-shared", "TX_SHARED", 24, 24, "TT")

    results = annotate_phased_variants(
        [
            make_isovar_result(
                v1,
                {"r12"},
                {"p12"},
                amino_acids="AA-1",
                cdna_sequences={"cdna-1"},
                transcript_ids={"tx-1", "tx-shared"},
                transcript_names={"TX1", "TX_SHARED"},
                known_somatic_transcript_edits=(somatic_edit_v1,),
                known_germline_transcript_edits=(germline_edit,),
                unexplained_transcript_edits=(unexplained_edit_v1,),
            ),
            make_isovar_result(
                v2,
                {"r12", "r23"},
                {"p12", "p23"},
                amino_acids="AA-2",
                cdna_sequences={"cdna-2"},
                transcript_ids={"tx-2", "tx-shared"},
                transcript_names={"TX2", "TX_SHARED"},
                known_somatic_transcript_edits=(somatic_edit_v2,),
                known_germline_transcript_edits=(germline_edit,),
                unexplained_transcript_edits=(unexplained_edit_v2,),
            ),
            make_isovar_result(
                v3,
                {"r23"},
                {"p23"},
                amino_acids="AA-3",
                cdna_sequences={"cdna-3"},
                transcript_ids={"tx-2"},
                transcript_names={"TX2"},
                known_somatic_transcript_edits=(somatic_edit_v3,),
            ),
        ],
        min_shared_fragments_for_phasing=1,
    )

    results_by_variant = {result.variant: result for result in results}
    expected_group_variants = (v1, v2, v3)
    expected_phase_read_names = frozenset({"r12", "r23"})
    expected_phase_protein_read_names = frozenset({"p12", "p23"})
    expected_cdna_sequences = ("cdna-1", "cdna-2", "cdna-3")
    expected_protein_sequences = ("AA-1", "AA-2", "AA-3")
    expected_transcript_ids = ("tx-1", "tx-2", "tx-shared")
    expected_transcript_names = ("TX1", "TX2", "TX_SHARED")
    expected_known_somatic_edits = (
        somatic_edit_v1,
        somatic_edit_v2,
        somatic_edit_v3,
    )
    expected_known_germline_edits = (germline_edit,)
    expected_unexplained_edits = (
        unexplained_edit_v1,
        unexplained_edit_v2,
    )
    expected_supporting_neighbors = {
        v1: {v2},
        v2: {v1, v3},
        v3: {v2},
    }

    for variant, result in results_by_variant.items():
        supporting_group = result.phase_group_from_supporting_reads
        protein_group = result.phase_group_from_protein_sequence

        assert supporting_group is not None
        assert protein_group is not None

        eq_(supporting_group.somatic_variants, expected_group_variants)
        eq_(supporting_group.germline_variants, ())
        eq_(supporting_group.supporting_read_names, expected_phase_read_names)
        eq_(supporting_group.cdna_sequences, ())
        eq_(supporting_group.mutant_protein_sequences, ())
        eq_(supporting_group.transcript_ids, ())
        eq_(supporting_group.transcript_names, ())
        eq_(supporting_group.known_somatic_transcript_edits, ())
        eq_(supporting_group.known_germline_transcript_edits, ())
        eq_(supporting_group.unexplained_transcript_edits, ())

        eq_(protein_group.somatic_variants, expected_group_variants)
        eq_(protein_group.germline_variants, ())
        eq_(protein_group.supporting_read_names, expected_phase_protein_read_names)
        eq_(protein_group.cdna_sequences, expected_cdna_sequences)
        eq_(protein_group.mutant_protein_sequences, expected_protein_sequences)
        eq_(protein_group.transcript_ids, expected_transcript_ids)
        eq_(protein_group.transcript_names, expected_transcript_names)
        eq_(
            protein_group.known_somatic_transcript_edits,
            expected_known_somatic_edits,
        )
        eq_(
            protein_group.known_germline_transcript_edits,
            expected_known_germline_edits,
        )
        eq_(
            protein_group.unexplained_transcript_edits,
            expected_unexplained_edits,
        )

        eq_(
            result.phased_variants_in_supporting_reads,
            expected_supporting_neighbors[variant],
        )
        eq_(
            result.phased_variants_in_protein_sequence,
            expected_supporting_neighbors[variant],
        )


def test_annotate_phased_variants_leaves_singletons_without_phase_group():
    variant = Variant("1", 20, "C", "A", normalize_contig_names=False)
    result = annotate_phased_variants(
        [make_isovar_result(variant, {"solo"}, {"solo-protein"})],
        min_shared_fragments_for_phasing=1,
    )[0]

    assert result.phase_group_from_supporting_reads is None
    assert result.phase_group_from_protein_sequence is None
    eq_(result.phased_variants_in_supporting_reads, set())
    eq_(result.phased_variants_in_protein_sequence, set())


def scoped_result(position, groups, segments=(64,), protein_groups=None):
    variant = Variant("1", position, "A", "C", normalize_contig_names=False)
    result = make_isovar_result(variant, (), {"shared"})
    reads = [AlleleRead(
        prefix="A", allele="C", suffix="T", name="shared",
        source_alignments=(((group, "shared", segment), (0, 0, "100M", False)),),
    ) for group in groups for segment in segments]
    result.read_evidence.alt_reads.extend(reads)
    result.top_protein_sequence.supporting_reads = {
        read for read in reads
        if protein_groups is None or read.source_alignments[0][0][0] in protein_groups
    }
    return result


@pytest.mark.parametrize("groups", [("a", "b"), ("", "a")])
def test_different_read_groups_cannot_phase_identical_names(groups):
    inputs = [scoped_result(10, [groups[0]]), scoped_result(20, [groups[1]])]
    # Public helper contracts still expose unmodified strings, not scoped IDs.
    for helper in (create_variant_to_alt_read_names_dict,
                   create_variant_to_protein_sequence_read_names_dict):
        assert list(helper(inputs).values()) == [{"shared"}, {"shared"}]
    for result in annotate_phased_variants(inputs, min_shared_fragments_for_phasing=1):
        assert result.alt_read_names == {"shared"}
        assert not result.phased_variants_in_supporting_reads
        assert not result.phased_variants_in_protein_sequence
        assert result.phase_group_from_supporting_reads is None
        assert result.phase_group_from_protein_sequence is None


@pytest.mark.parametrize("group", ["a", ""])
def test_complementary_mates_in_same_group_still_phase(group):
    inputs = [scoped_result(10, [group], (64,)), scoped_result(20, [group], (128,))]
    for result in annotate_phased_variants(inputs, min_shared_fragments_for_phasing=1):
        assert result.phased_variants_in_supporting_reads == {
            r.variant for r in inputs if r.variant != result.variant}
        assert result.phased_variants_in_protein_sequence == result.phased_variants_in_supporting_reads
        for phase_group in (result.phase_group_from_supporting_reads,
                            result.phase_group_from_protein_sequence):
            assert phase_group.supporting_read_names == frozenset({"shared"})


def test_read_groups_count_separately_but_mates_and_duplicate_placements_do_not():
    for groups, expected in [(["a"], False), (["a", "b"], True)]:
        inputs = [scoped_result(pos, groups, (64, 128, 64)) for pos in (10, 20)]
        for result in annotate_phased_variants(inputs, min_shared_fragments_for_phasing=2):
            assert bool(result.phased_variants_in_supporting_reads) == expected
            assert bool(result.phased_variants_in_protein_sequence) == expected
            if expected:
                assert result.phase_group_from_protein_sequence.supporting_read_names == {"shared"}


def test_protein_phasing_uses_its_actual_subset_not_matching_alt_names():
    inputs = [scoped_result(10, ["a", "b"], protein_groups=["a"]),
              scoped_result(20, ["a", "b"], protein_groups=["b"])]
    for result in annotate_phased_variants(inputs, min_shared_fragments_for_phasing=1):
        assert result.phased_variants_in_supporting_reads
        assert not result.phased_variants_in_protein_sequence
        assert result.phase_group_from_protein_sequence is None


def test_missing_metadata_does_not_match_a_collected_unscoped_read():
    legacy = make_isovar_result(Variant("1", 10, "A", "C"), {"shared"}, {"shared"})
    collected = scoped_result(20, [""])
    for result in annotate_phased_variants([legacy, collected], min_shared_fragments_for_phasing=1):
        assert not result.phased_variants_in_supporting_reads
        assert not result.phased_variants_in_protein_sequence


def test_public_phase_group_helper_still_accepts_plain_names():
    variants = [Variant("1", pos, "A", "C") for pos in (10, 20)]
    groups = create_phase_groups({v: {"shared"} for v in variants}, 1)
    assert set(groups) == set(variants)
    assert all(g.supporting_read_names == {"shared"} for g in groups.values())


def test_read_collector_and_adapter_do_not_link_different_groups():
    from isovar.read_collector import ReadCollector
    from isovar.read_phasing import IsovarReadPhasing
    from tests.mock_objects import MockAlignmentFile, make_pysam_read

    records = []
    for group, sequence in [("a", "ACCGTGATCG"), ("b", "ACCTTGAACG")]:
        read = make_pysam_read(sequence, "10M", name="shared", reference_start=0)
        read.set_tag("RG", group)
        records.append(read)
    bam = MockAlignmentFile(["1"], records)
    variants = [Variant("1", 4, "T", "G"), Variant("1", 8, "T", "A")]
    inputs = [IsovarResult(v, ReadCollector().read_evidence_for_variant(v, bam), None)
              for v in variants]
    results = annotate_phased_variants(inputs, min_shared_fragments_for_phasing=1)
    adapter = IsovarReadPhasing(results)
    assert all(r.num_alt_fragments == 1 and r.alt_read_names == {"shared"} for r in results)
    assert all(adapter.has_evidence(v) and not adapter.partners_in_cis(v) for v in variants)
