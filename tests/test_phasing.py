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

from varcode import Variant
from varcode.mutant_transcript import TranscriptEdit

from isovar.allele_read import AlleleRead
from isovar.isovar_result import IsovarResult
from isovar.phasing import annotate_phased_variants
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
