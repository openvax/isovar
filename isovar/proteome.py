"""Adapt Isovar's small-variant protein hypotheses to ORFdata Proteome.

The adapter copies recorded RNA observations and translations; it does not
infer protein expression, complete partial sequences, or select a winner.
"""

from collections.abc import Mapping

from orfdata import ProducerInfo, ProducerLink, ProducerRecord, Proteome, ProteomeEntry
from orfdata import ReferenceContext as ProteomeReferenceContext

from .protein_hypotheses import SCHEMA, export_protein_hypotheses
from .rna_evidence import content_identifier


def _record(kind, native_id, payload, sample_id, source, **kwargs):
    return ProducerRecord(
        namespace="isovar", schema_version="v2", kind=kind, native_id=native_id,
        exporter="isovar", sample_id=sample_id, run_id=source, payload=payload, **kwargs)


def _support_link(support, target, evidence_records):
    if not isinstance(support, Mapping):
        raise ValueError("RNA support must be a JSON object")
    evidence_id = support.get("evidence_set_id")
    if evidence_id is None:
        return None
    evidence = evidence_records.get(evidence_id)
    if evidence is None:
        raise ValueError(f"Unknown RNA evidence set: {evidence_id}")
    return ProducerLink(evidence.key, target.key, "record", "supports_rna_sequence")


def _event_status(proteins, complete):
    if complete is False:
        return "capped"
    if not proteins and complete is True:
        return "zero_product"
    return "retained" if complete is True else "completeness_unknown"


def from_protein_hypotheses(export, *, reference):
    """Convert a portable Isovar v2 export without reads or a live genome.

    Parameters
    ----------
    export : dict
        Output of :func:`isovar.export_protein_hypotheses`, or its parsed JSON.
    reference : orfdata.ReferenceContext
        Explicit taxon, build and annotation asserted by the caller. The v2
        export only records a build name, so it cannot verify the annotation
        release. A contradictory recorded build is rejected.

    Returns
    -------
    orfdata.Proteome
        One entry per retained protein hypothesis, with separate records for
        events, hypotheses, translations and scoped RNA evidence. A no-protein
        event still has an event record. RNA support is not protein evidence.
    """
    if not isinstance(export, Mapping) or export.get("schema") != SCHEMA:
        raise ValueError(f"Expected {SCHEMA} export")
    if not isinstance(reference, ProteomeReferenceContext) or not reference.is_genomic:
        raise ValueError("A genomic ORFdata ReferenceContext is required")
    sample_id, source = export.get("sample_id"), export.get("source")
    if any(not isinstance(value, str) or not value.strip() for value in (sample_id, source)):
        raise ValueError("Isovar sample_id and source must be nonempty strings")
    if export.get("evidence_scope") != [sample_id, source]:
        raise ValueError("Isovar evidence scope does not match sample_id and source")
    if export.get("interval_convention") != "zero_based_half_open":
        raise ValueError("Unsupported Isovar interval convention")
    events, evidence_sets = export.get("events"), export.get("evidence_sets")
    if not isinstance(events, list) or not isinstance(evidence_sets, Mapping):
        raise ValueError("Isovar events and evidence_sets are required")
    version = export.get("isovar_version")
    if not isinstance(version, str) or not version:
        raise ValueError("Isovar version is required")

    metadata = {key: value for key, value in export.items() if key not in ("events", "evidence_sets")}
    exported = _record("export", content_identifier("export", export), metadata, sample_id, source)
    records = [exported]
    links = []
    entries = []
    evidence_records = {}
    for evidence_id, payload in sorted(evidence_sets.items()):
        if not isinstance(evidence_id, str) or not isinstance(payload, Mapping):
            raise ValueError("RNA evidence sets need string IDs and JSON objects")
        record = _record(
            "evidence", evidence_id, payload, sample_id, source, support="rna_sequence",
            scope={"sample": sample_id, "read_set": source})
        evidence_records[evidence_id] = record
        records.append(record)

    for event_index, event in enumerate(events):
        if not isinstance(event, Mapping) or not isinstance(event.get("protein_hypotheses"), list):
            raise ValueError("Each Isovar event needs protein_hypotheses")
        event_id = event.get("event_id")
        if not isinstance(event_id, str) or not event_id:
            raise ValueError("Each Isovar event needs an event_id")
        variant = event.get("variant")
        if not isinstance(variant, Mapping):
            raise ValueError("Each Isovar event needs a variant")
        declared_build = variant.get("reference_name")
        if declared_build is not None and declared_build != reference.genome_build:
            raise ValueError(f"Isovar build {declared_build!r} conflicts with {reference.genome_build!r}")
        proteins = event["protein_hypotheses"]
        complete = event.get("protein_hypotheses_complete")
        if complete is not None and type(complete) is not bool:
            raise ValueError("protein_hypotheses_complete must be true, false or null")
        event_payload = {key: value for key, value in event.items() if key != "protein_hypotheses"}
        event_record = _record(
            "event", f"{event_id}:{event_index}", event_payload, sample_id, source,
            status=_event_status(proteins, complete), complete=complete,
            limits={"max_protein_sequences_per_variant": event.get("protein_sequence_limit")})
        records.append(event_record)
        links.append(ProducerLink(exported.key, event_record.key, "record", "contains_event"))
        allele_support = event.get("allele_support", {})
        if not isinstance(allele_support, Mapping):
            raise ValueError("allele_support must be a JSON object")
        if "alt" in allele_support:
            link = _support_link(allele_support["alt"], event_record, evidence_records)
            if link is not None:
                links.append(link)

        for protein_index, protein in enumerate(proteins):
            if not isinstance(protein, Mapping) or not isinstance(protein.get("translations"), list):
                raise ValueError("Each protein hypothesis needs translations")
            hypothesis_id, sequence = protein.get("hypothesis_id"), protein.get("amino_acids")
            if not isinstance(hypothesis_id, str) or not hypothesis_id:
                raise ValueError("Each protein hypothesis needs a hypothesis_id")
            native_id = f"{hypothesis_id}:event{event_index}:rank{protein_index + 1}"
            hypothesis_payload = {key: value for key, value in protein.items() if key != "translations"}
            hypothesis_record = _record(
                "hypothesis", native_id, hypothesis_payload, sample_id, source)
            records.append(hypothesis_record)
            links.append(ProducerLink(event_record.key, hypothesis_record.key, "record", "has_hypothesis"))
            n_terminus = protein.get("n_terminus")
            c_terminus = protein.get("c_terminus")
            if n_terminus not in ("annotated_start_codon", "local_window", "mixed"):
                raise ValueError(f"Unknown Isovar N terminus: {n_terminus!r}")
            if c_terminus not in ("stop_codon", "no_stop_in_context"):
                raise ValueError(f"Unknown Isovar C terminus: {c_terminus!r}")
            producer = ProducerInfo(
                origin="reference_cds", method="rna_assembly", tool="isovar",
                native_id=native_id, tool_version=version, exporter="isovar",
                sample_id=sample_id, run_id=source,
                gene_ids=tuple(protein.get("gene_ids", ())),
                gene_names=tuple(protein.get("gene_names", ())),
                transcript_ids=tuple(protein.get("transcript_ids", ())),
                reference=reference,
                start_complete=(True if n_terminus == "annotated_start_codon" else
                                False if n_terminus == "local_window" else None),
                stop_complete=c_terminus == "stop_codon")
            entry = ProteomeEntry.from_producer(sequence, taxon_id=reference.taxon_id, producer=producer)
            entries.append(entry)
            links.append(ProducerLink(hypothesis_record.key, entry.accession, "entry", "has_sequence"))
            link = _support_link(protein.get("rna_support"), hypothesis_record, evidence_records)
            if link is not None:
                links.append(link)
            for translation_index, translation in enumerate(protein["translations"]):
                if not isinstance(translation, Mapping):
                    raise ValueError("Each Isovar translation must be a JSON object")
                translation_id = translation.get("translation_id")
                if not isinstance(translation_id, str) or not translation_id:
                    raise ValueError("Each Isovar translation needs a translation_id")
                translation_record = _record(
                    "translation", f"{translation_id}:event{event_index}:protein{protein_index}:"
                    f"translation{translation_index}", translation, sample_id, source)
                records.append(translation_record)
                links.append(ProducerLink(hypothesis_record.key, translation_record.key,
                                          "record", "has_translation"))
                link = _support_link(translation.get("rna_support"), translation_record, evidence_records)
                if link is not None:
                    links.append(link)

    return Proteome.from_entries(entries, reference=reference, records=records, links=links)


def to_proteome(isovar_results, *, sample_id, source, reference, **export_options):
    """Export native Isovar results and adapt them using the offline mapping."""
    export = export_protein_hypotheses(
        isovar_results, sample_id=sample_id, source=source, **export_options)
    return from_protein_hypotheses(export, reference=reference)
