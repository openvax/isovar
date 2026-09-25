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

__version__ = "1.36.0"


from .allele_interpretations import reconcile_allele_interpretations
from .allele_read import AlleleRead
from .cell_evidence import cell_umi_allele_evidence
from .dataframe_helpers import isovar_results_to_dataframe
from .isovar_result import IsovarResult
from .locus_read import LocusRead
from .main import run_isovar
from .mutant_transcript_source import IsovarMutantTranscript
from .phase_group import PhaseGroup
from .protein_hypotheses import export_protein_hypotheses, read_group_metadata, write_protein_hypotheses
from .protein_sequence import ProteinSequence
from .protein_sequence_creator import ProteinSequenceCreator
from .read_collector import ReadCollector
from .read_end_inference import (
    Adapter, ReadEndProfile, ReadEndAnnotation, ReadSequenceView,
    infer_read_ends, read_sequence_view_from_alignment,
)
from .read_evidence import ReadEvidence
from .read_phasing import IsovarReadPhasing
from .rna_evidence import union_rna_support
from .transcript_assembly_edit import TranscriptAssemblyEdit
from .variant_orf import VariantORF
from .variant_sequence import VariantSequence
from .variant_sequence_creator import VariantSequenceCreator
from .fusion import FusionBlock, FusionBreakpoint, FusionRead, FusionReference, FusionTranscript, reconstruct_fusion
from .sv_rna import reconstruct_sv_rna
from .sv_rna_orf_export import export_sv_rna_orfs, write_sv_rna_orfs
from .sv_rna_comparison import compare_sv_rna_predictions
from .orf_inclusion import annotate_orf_inclusion
from .orf_start import annotate_orf_start, summarize_orf_start_evidence


__all__ = [
    "run_isovar",
    "reconcile_allele_interpretations",
    "cell_umi_allele_evidence",
    "export_protein_hypotheses",
    "write_protein_hypotheses",
    "read_group_metadata",
    "union_rna_support",
    "FusionBlock",
    "FusionBreakpoint",
    "FusionRead",
    "FusionReference",
    "FusionTranscript",
    "reconstruct_fusion",
    "reconstruct_sv_rna",
    "export_sv_rna_orfs",
    "write_sv_rna_orfs",
    "compare_sv_rna_predictions",
    "annotate_orf_inclusion",
    "annotate_orf_start",
    "summarize_orf_start_evidence",
    "isovar_results_to_dataframe",
    "AlleleRead",
    "PhaseGroup",
    "IsovarMutantTranscript",
    "IsovarReadPhasing",
    "IsovarResult",
    "LocusRead",
    "ProteinSequence",
    "ProteinSequenceCreator",
    "ReadCollector",
    "Adapter",
    "ReadEndProfile",
    "ReadEndAnnotation",
    "ReadSequenceView",
    "infer_read_ends",
    "read_sequence_view_from_alignment",
    "ReadEvidence",
    "TranscriptAssemblyEdit",
    "VariantORF",
    "VariantSequence",
    "VariantSequenceCreator",
]
