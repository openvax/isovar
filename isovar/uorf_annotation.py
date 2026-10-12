"""Compatibility wrapper for the human upstream-ORF catalogue."""

from .orf_annotation import annotate_sv_orf_evidence
from .uorf_catalog import load_uorf_catalog


def annotate_sv_uorf_evidence(result, transcripts, *, catalog=None):
    annotated = annotate_sv_orf_evidence(
        result, transcripts, catalog=load_uorf_catalog() if catalog is None else catalog)
    annotated["uorf_catalogue_metadata"] = annotated.pop("orf_catalogue_metadata")
    for path in annotated.get("paths", ()):
        for candidate in path.get("exploratory_orfs", {}).get("candidates", ()):
            candidate["reference_uorf_annotations"] = candidate.pop("reference_orf_annotations")
            candidate["reference_uorf_lookup"] = candidate.pop("reference_orf_lookup")
    return annotated
