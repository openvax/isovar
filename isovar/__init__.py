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

from importlib import import_module

__version__ = "1.95.0"

# Lazily resolved top-level public API so importing ``isovar`` (or its CLI)
# stays fast and avoids loading heavy runtime dependencies until needed.
_LAZY_EXPORTS = {
    "reconcile_allele_interpretations": (".allele_interpretations", "reconcile_allele_interpretations"),
    "AlleleRead": (".allele_read", "AlleleRead"),
    "BaseQualityPolicy": (".base_quality", "BaseQualityPolicy"),
    "AlleleQualityFootprint": (".base_quality", "AlleleQualityFootprint"),
    "cell_umi_allele_evidence": (".cell_evidence", "cell_umi_allele_evidence"),
    "UnresolvedContigError": (".contig_resolution", "UnresolvedContigError"),
    "CallerCoordinateClaim": (".sv_sources", "CallerCoordinateClaim"),
    "SVSourceRecord": (".sv_sources", "SVSourceRecord"),
    "read_sv_source": (".sv_sources", "read_sv_source"),
    "sv_source_report": (".sv_sources", "sv_source_report"),
    "PreparedSVSourceRecord": (".sv_source_preparation", "PreparedSVSourceRecord"),
    "SVSourceView": (".sv_source_preparation", "SVSourceView"),
    "prepare_sv_source_record": (".sv_source_preparation", "prepare_sv_source_record"),
    "isovar_results_to_dataframe": (".dataframe_helpers", "isovar_results_to_dataframe"),
    "IsovarResult": (".isovar_result", "IsovarResult"),
    "LocusRead": (".locus_read", "LocusRead"),
    "run_isovar": (".main", "run_isovar"),
    "IsovarMutantTranscript": (".mutant_transcript_source", "IsovarMutantTranscript"),
    "PhaseGroup": (".phase_group", "PhaseGroup"),
    "export_protein_hypotheses": (".protein_hypotheses", "export_protein_hypotheses"),
    "read_group_metadata": (".protein_hypotheses", "read_group_metadata"),
    "write_protein_hypotheses": (".protein_hypotheses", "write_protein_hypotheses"),
    "from_protein_hypotheses": (".proteome", "from_protein_hypotheses"),
    "to_proteome": (".proteome", "to_proteome"),
    "from_rna_result": (".rna_proteome", "from_rna_result"),
    "export_patient_background_pairs": (".patient_background", "export_patient_background_pairs"),
    "reconstruct_cell_groups": (".cell_reconstruction", "reconstruct_cell_groups"),
    "ProteinSequence": (".protein_sequence", "ProteinSequence"),
    "ProteinSequenceCreator": (".protein_sequence_creator", "ProteinSequenceCreator"),
    "ReadCollector": (".read_collector", "ReadCollector"),
    "Adapter": (".read_end_inference", "Adapter"),
    "ReadEndProfile": (".read_end_inference", "ReadEndProfile"),
    "ReadEndAnnotation": (".read_end_inference", "ReadEndAnnotation"),
    "ReadSequenceView": (".read_end_inference", "ReadSequenceView"),
    "ReadEndInferenceCache": (".read_end_inference", "ReadEndInferenceCache"),
    "infer_read_ends": (".read_end_inference", "infer_read_ends"),
    "read_sequence_view_from_alignment": (".read_end_inference", "read_sequence_view_from_alignment"),
    "ReadLibraryProvenance": (".read_end_profiles", "ReadLibraryProvenance"),
    "ReadProcessingStep": (".read_end_profiles", "ReadProcessingStep"),
    "read_end_profile_from_fasta": (".read_end_profiles", "read_end_profile_from_fasta"),
    "ReadEvidence": (".read_evidence", "ReadEvidence"),
    "GenomicQueryEvidence": (".read_sequence_origin", "GenomicQueryEvidence"),
    "genomic_query_evidence_from_alignment": (".read_sequence_origin", "genomic_query_evidence_from_alignment"),
    "PhasingCalibration": (".phasing_calibration", "PhasingCalibration"),
    "PhasingControl": (".phasing_calibration", "PhasingControl"),
    "calibrate_phasing_error_rates": (".phasing_calibration", "calibrate_phasing_error_rates"),
    "IsovarReadPhasing": (".read_phasing", "IsovarReadPhasing"),
    "union_rna_support": (".rna_evidence", "union_rna_support"),
    "RNAFeature": (".rna_feature_support", "RNAFeature"),
    "RNASamplingGeometry": (".rna_feature_support", "RNASamplingGeometry"),
    "RNASamplingProfile": (".rna_feature_support", "RNASamplingProfile"),
    "effective_rna_feature_lengths": (".rna_feature_support", "effective_rna_feature_lengths"),
    "measure_rna_feature": (".rna_feature_support", "measure_rna_feature"),
    "validate_rna_feature_support": (".rna_feature_support", "validate_rna_feature_support"),
    "normalize_rna_feature_support": (".rna_feature_support", "normalize_rna_feature_support"),
    "measure_rna_path_features": (".rna_feature_measurement", "measure_rna_path_features"),
    "RNAOriginHypothesis": (".rna_origin", "RNAOriginHypothesis"),
    "RNAOriginPolicy": (".rna_origin", "RNAOriginPolicy"),
    "assess_rna_origin": (".rna_origin", "assess_rna_origin"),
    "TranscriptAssemblyEdit": (".transcript_assembly_edit", "TranscriptAssemblyEdit"),
    "VariantORF": (".variant_orf", "VariantORF"),
    "VariantRNAAssembly": (".variant_read_assessment", "VariantRNAAssembly"),
    "VariantReadAssessment": (".variant_read_assessment", "VariantReadAssessment"),
    "VariantSequence": (".variant_sequence", "VariantSequence"),
    "VariantSequenceCreator": (".variant_sequence_creator", "VariantSequenceCreator"),
    "FUSION_REASON_ALIASES": (".fusion", "FUSION_REASON_ALIASES"),
    "FusionBlock": (".fusion", "FusionBlock"),
    "FusionBreakpoint": (".fusion", "FusionBreakpoint"),
    "FusionRead": (".fusion", "FusionRead"),
    "FusionReference": (".fusion", "FusionReference"),
    "FusionTranscript": (".fusion", "FusionTranscript"),
    "normalize_fusion_reason": (".fusion", "normalize_fusion_reason"),
    "reconstruct_fusion": (".fusion", "reconstruct_fusion"),
    "reconstruct_sv_rna": (".sv_rna", "reconstruct_sv_rna"),
    "export_sv_rna_orfs": (".sv_rna_orf_export", "export_sv_rna_orfs"),
    "write_sv_rna_orfs": (".sv_rna_orf_export", "write_sv_rna_orfs"),
    "compare_sv_rna_predictions": (".sv_rna_comparison", "compare_sv_rna_predictions"),
    "assess_sv_rna_paths": (".sv_rna_assessment", "assess_sv_rna_paths"),
    "hypotheses_from_varcode": (".sv_rna_reconciliation", "hypotheses_from_varcode"),
    "reconcile_sv_rna_hypotheses": (".sv_rna_reconciliation", "reconcile_sv_rna_hypotheses"),
    "run_sv_rna_audit": (".sv_rna_audit", "run_sv_rna_audit"),
    "annotate_orf_inclusion": (".orf_inclusion", "annotate_orf_inclusion"),
    "annotate_orf_start": (".orf_start", "annotate_orf_start"),
    "summarize_orf_start_evidence": (".orf_start", "summarize_orf_start_evidence"),
    "UORFCatalog": (".uorf_catalog", "UORFCatalog"),
    "UORFRecord": (".uorf_catalog", "UORFRecord"),
    "UORFPeptide": (".uorf_catalog", "UORFPeptide"),
    "load_uorf_catalog": (".uorf_catalog", "load_uorf_catalog"),
    "InitiationContext": (".initiation_context", "InitiationContext"),
    "annotate_initiation_context": (".initiation_context", "annotate_initiation_context"),
    "IsoformTarget": (".antigen_reference", "IsoformTarget"),
    "gdc_genome": (".antigen_reference", "gdc_genome"),
    "gdc_viral_panel": (".antigen_reference", "gdc_viral_panel"),
    "isoform_targets_from_genome": (".antigen_reference", "isoform_targets_from_genome"),
    "validate_alignment_reference": (".antigen_reference", "validate_alignment_reference"),
    "reconstruct_antigen_isoforms": (".antigen_isoforms", "reconstruct_antigen_isoforms"),
    "AntigenReference": (".antigen_targets", "AntigenReference"),
    "antigen_gene_selection": (".antigen_targets", "antigen_gene_selection"),
    "prepare_antigen_reference": (".antigen_targets", "prepare_antigen_reference"),
    "reconstruct_viral_sequences": (".viral_sequences", "reconstruct_viral_sequences"),
    "write_viral_fasta": (".viral_sequences", "write_viral_fasta"),
    "write_viral_proteins": (".viral_coding", "write_viral_proteins"),
    "quantify_antigen_isoforms": (".isoform_quantification", "quantify_antigen_isoforms"),
    "add_isoform_abundances": (".isoform_quantification", "add_isoform_abundances"),
    "write_isoform_fastas": (".isoform_quantification", "write_isoform_fastas"),
    "annotate_sv_uorf_evidence": (".uorf_annotation", "annotate_sv_uorf_evidence"),
    "annotate_sv_orf_evidence": (".orf_annotation", "annotate_sv_orf_evidence"),
    "reconstruct_variant_orfs": (".variant_rna_orfs", "reconstruct_variant_orfs"),
    "variant_orf_translation_table": (".variant_rna_orfs", "variant_orf_translation_table"),
    "prepare_rna_library": (".rna_library", "prepare_rna_library"),
    "PreparedRnaLibrary": (".rna_library", "PreparedRnaLibrary"),
    "RnaLibraryIndex": (".rna_library", "RnaLibraryIndex"),
    "PathEvidenceQuery": (".sv_path_routing", "PathEvidenceQuery"),
    "PathEvidenceIndex": (".sv_path_routing", "PathEvidenceIndex"),
    "route_path_evidence": (".sv_path_routing", "route_path_evidence"),
    "NormalEvidencePolicy": (".normal_evidence", "NormalEvidencePolicy"),
    "NormalSource": (".normal_sources", "NormalSource"),
    "assess_matched_normal": (".normal_sources", "assess_matched_normal"),
    "run_normal_evidence_manifest": (".normal_sources", "run_normal_evidence_manifest"),
    "normal_background_regions": (".normal_regions", "normal_background_regions"),
    "reconstruct_variant_rna_paths": (".variant_rna_paths", "reconstruct_variant_rna_paths"),
    "variant_rna_path_input_from_dict": (".variant_rna_paths", "variant_rna_path_input_from_dict"),
}

__all__ = list(_LAZY_EXPORTS)


def __getattr__(name):
    if name in _LAZY_EXPORTS:
        module_name, attr_name = _LAZY_EXPORTS[name]
        value = getattr(import_module(module_name, __name__), attr_name)
        globals()[name] = value
        return value
    raise AttributeError("module %r has no attribute %r" % (__name__, name))


def __dir__():
    return sorted(set(globals()) | set(__all__))
