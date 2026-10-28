"""Preserve RNA producer hypotheses and evidence in ORFdata snapshots."""

from collections.abc import Mapping
from hashlib import sha256

from orfdata import ProducerInfo, ProducerLink, Proteome, ReferenceContext

from .genetic_code import standard_genetic_code, translate_cdna
from .proteome import (
    _PROTEIN_TARGET, _Protein, _record, _support_link, _validate_evidence,
    from_protein_hypotheses,
)
from .rna_evidence import canonical_json
from .sv_rna_orf_export import export_sv_rna_orfs


def _rows(value, name):
    if not isinstance(value, list) or any(not isinstance(row, Mapping) for row in value):
        raise ValueError(f"{name} must be an array of JSON objects")
    return value


def _indexed(rows, field):
    result = {}
    for row in rows:
        key = row.get(field)
        if not isinstance(key, str) or not key or key in result:
            raise ValueError(f"Missing or duplicate {field}")
        result[key] = row
    return result


def _translation(sequence, protein, stop, *, table=1, start=False):
    if not isinstance(sequence, str) or not isinstance(protein, str) or not protein or type(stop) is not bool:
        raise ValueError("A protein candidate needs nucleotide/protein sequences and stop status")
    if len(sequence) != 3 * (len(protein) + int(stop)):
        raise ValueError("RNA coding interval and protein length disagree")
    amino, observed_stop = translate_cdna(sequence, translation_table=table)
    if start and amino:
        amino = 'M' + amino[1:]
    if amino != protein or observed_stop != stop:
        raise ValueError("RNA coding sequence and protein translation disagree")


def _coding_sequence(sequence, interval):
    if (not isinstance(sequence, str) or not isinstance(interval, (list, tuple)) or len(interval) != 2
            or any(type(i) is not int for i in interval) or not 0 <= interval[0] < interval[1] <= len(sequence)):
        raise ValueError("Invalid RNA coding interval")
    return sequence[interval[0]:interval[1]]


class _Snapshot:
    """Build one source-qualified record graph without selecting a protein."""

    def __init__(self, value, reference, source):
        self.reference, self.sample, self.source = reference, value.get('sample_id'), source
        if any(not isinstance(v, str) or not v.strip() for v in (self.sample, source)):
            raise ValueError("Isovar sample_id and source must be nonempty strings")
        self.schema = value['schema'].removeprefix('isovar.')
        self.digest = sha256(canonical_json(value).encode()).hexdigest()
        self.version = value.get('producer', {}).get('version') or value.get('isovar_version')
        self.records, self.links, self.proteins, self.evidence = [], [], [], {}

    def record(self, kind, native_id, payload, **kwargs):
        record = _record(kind, self.digest + ':' + native_id, payload, self.sample, self.source,
                         schema_version=self.schema, **kwargs)
        self.records.append(record)
        return record

    def link(self, parent, child, relation, **payload):
        self.links.append(ProducerLink(parent.key, child.key, 'record', relation, payload))

    def add_evidence(self, evidence_id, payload):
        _validate_evidence(evidence_id, payload, self.sample, self.source)
        prior = self.evidence.get(evidence_id)
        if prior is not None:
            if prior.to_dict()['payload'] != dict(payload):
                raise ValueError("Conflicting RNA evidence records with the same identity")
            return prior
        record = self.record('evidence', evidence_id, payload, support='rna_sequence',
                             scope={'sample': self.sample, 'read_set': self.source})
        self.evidence[evidence_id] = record
        return record

    def support(self, support, target, *, inline=False):
        if inline and isinstance(support, Mapping) and support.get('evidence_set_id') is not None:
            keys = ('evidence_set_id', 'evidence_scope', 'reads', 'fragments', 'read_ids', 'fragment_ids')
            self.add_evidence(support['evidence_set_id'], {k: support.get(k) for k in keys})
        link = _support_link(support, target, self.evidence)
        if link is not None:
            self.links.append(link)

    def protein(self, record, native_id, amino, origin, *, start, stop, **identities):
        producer = ProducerInfo(
            origin=origin, method='rna_assembly', tool='isovar', native_id=native_id,
            tool_version=self.version, exporter='isovar', sample_id=self.sample, run_id=self.source,
            reference=self.reference, source_sha256=self.digest,
            start_complete=start, stop_complete=stop, **identities)
        protein = _Protein.from_producer(amino, taxon_id=self.reference.taxon_id, producer=producer)
        self.proteins.append(protein)
        self.links.append(ProducerLink(record.key, protein.accession, _PROTEIN_TARGET, 'has_sequence'))

    def finish(self):
        return Proteome.from_proteins(self.proteins, reference=self.reference,
                                     records=self.records, links=self.links)


def _sv(value, reference):
    if value.get('interval_convention') != 'zero_based_half_open':
        raise ValueError("Unsupported Isovar interval convention")
    candidates = _rows(value.get('candidates'), 'SV candidates')
    _indexed(candidates, 'candidate_id')
    graph = _Snapshot(value, reference, value.get('source'))
    event = graph.record('event', value['event_id'], {k: v for k, v in value.items() if k != 'candidates'},
                         status='limited' if value['limitations'] else 'retained' if candidates else 'zero_product')
    for index, candidate in enumerate(candidates):
        cid = candidate['candidate_id']
        dna, amino, stop = candidate['nucleotide_sequence'], candidate['amino_acids'], candidate['ends_with_stop_codon']
        if not dna.startswith('ATG'):
            raise ValueError("Exploratory SV ORFs must retain their ATG start")
        _translation(dna, amino, stop, start=True)
        occurrences = _rows(candidate.get('occurrences'), 'SV occurrences')
        if not occurrences or candidate.get('initiation_observed') is not False or candidate.get('translation_observed') is not False:
            raise ValueError("SV RNA candidates require occurrences and unobserved initiation/translation")
        record = graph.record('hypothesis', cid, {k: v for k, v in candidate.items() if k != 'occurrences'},
                              status='exploratory')
        graph.link(event, record, 'has_hypothesis', index=index)
        graph.support(candidate['rna_support'], record, inline=True)
        transcripts = set()
        measured = []
        for i, occurrence in enumerate(occurrences):
            lo, hi = occurrence['query_interval']
            if type(lo) is not int or type(hi) is not int or lo < 0 or hi - lo != len(dna):
                raise ValueError("SV occurrence interval disagrees with its nucleotide sequence")
            placed = graph.record('occurrence', cid + ':' + str(i), occurrence,
                                  status=occurrence['frame_status'])
            graph.link(record, placed, 'has_occurrence', index=i)
            if 'rna_support' in occurrence:
                graph.support(occurrence['rna_support'], placed, inline=True)
                scoped = occurrence['rna_support']
                if scoped.get('scope') != 'full_interval_junction_linked_witnesses_for_occurrence':
                    raise ValueError("SV occurrence RNA support must be measured for that occurrence")
                if bool(occurrence['witnesses']) != bool(scoped['reads']):
                    raise ValueError("SV occurrence witnesses and RNA support disagree")
                for ids in ('read_ids', 'fragment_ids'):
                    if not set(scoped[ids]) <= set(candidate['rna_support'][ids]):
                        raise ValueError("SV occurrence support lies outside the candidate union")
                measured.append(scoped)
            transcripts.update(c['transcript_id'] for c in occurrence['reference_comparisons'] if c.get('transcript_id'))
        if len(measured) == len(occurrences):
            for ids in ('read_ids', 'fragment_ids'):
                if set(candidate['rna_support'][ids]) != {value for support in measured for value in support[ids]}:
                    raise ValueError("SV candidate union disagrees with occurrence support")
        graph.protein(record, cid, amino, 'sv_orf',
                      start=True if all(o['annotated_start'] for o in occurrences) else None,
                      stop=stop, transcript_ids=tuple(sorted(transcripts)))
    return graph.finish()


def _variant(value, reference):
    candidates = _rows(value.get('protein_candidates'), 'Protein candidates')
    _indexed(candidates, 'candidate_id')
    models = _indexed(_rows(value.get('models'), 'ORF models'), 'model_id')
    rna = _indexed(_rows(value.get('rna_candidates'), 'RNA candidates'), 'rna_candidate_id')
    graph = _Snapshot(value, reference, value.get('source'))
    excluded = {'protein_candidates', 'models', 'rna_candidates', 'evidence_sets'}
    event = graph.record('event', 'variant', {k: v for k, v in value.items() if k not in excluded},
                         status=value['assembly_state'])
    if not isinstance(value.get('evidence_sets'), Mapping):
        raise ValueError("Variant RNA ORFs require evidence_sets")
    for eid, evidence in value['evidence_sets'].items():
        graph.add_evidence(eid, evidence)
    for key in ('allele_support', 'reference_incompatible_support'):
        graph.support(value[key], event)
    model_records, rna_records = {}, {}
    for i, (mid, model) in enumerate(models.items()):
        context = ReferenceContext(**model['reference']['reference_context'])
        if context.key != reference.key:
            raise ValueError("ORF model reference context conflicts with the supplied reference")
        record = graph.record('model', mid, model, status=model['state'])
        model_records[mid] = record
        graph.link(event, record, 'has_model', index=i)
    for i, (rid, row) in enumerate(rna.items()):
        record = graph.record('rna_candidate', rid, row, status=row['state'])
        rna_records[rid] = record
        graph.link(event, record, 'has_rna_candidate', index=i)
        graph.support(row['rna_support'], record)
    for index, candidate in enumerate(candidates):
        cid, mid, rid = candidate['candidate_id'], candidate['model_id'], candidate['rna_candidate_id']
        if mid not in models or rid not in rna or cid not in models[mid]['candidate_ids']:
            raise ValueError("Protein candidate has an unknown model/RNA linkage")
        model = models[mid]['reference']
        amino, stop, start = candidate['protein_sequence'], candidate['ends_with_stop_codon'], candidate['starts_at_reference_start']
        if type(start) is not bool or candidate.get('initiation_observed') is not False or candidate.get('translation_observed') is not False:
            raise ValueError("Variant RNA candidates must retain unobserved initiation/translation and start status")
        dna = _coding_sequence(candidate['nucleotide_sequence'], candidate['coding_interval'])
        if start and dna[:3] != model['start_codon']:
            raise ValueError("RNA start codon disagrees with the selected ORF model")
        _translation(dna, amino, stop, table=model['translation_table'], start=start)
        record = graph.record('hypothesis', cid, candidate, status=candidate['evidence_state'])
        graph.link(event, record, 'has_hypothesis', index=index)
        graph.link(model_records[mid], record, 'supplies_frame')
        graph.link(rna_records[rid], record, 'has_translation')
        graph.support(candidate['rna_support'], record)
        for key in ('single_read_span_support', 'allele_observation_span_support'):
            if key in candidate['coding_support']:
                graph.support(candidate['coding_support'][key], record)
        gene = model.get('source_gene_id') if model.get('on_source_transcript') else None
        if model['kind'] not in ('annotated', 'extra_orf'):
            raise ValueError("Unknown ORF model kind")
        graph.protein(record, cid, amino, 'reference_cds' if model['kind'] == 'annotated' else 'nc_orf',
                      start=start, stop=stop, gene_ids=(gene,) if gene else (),
                      gene_names=(candidate['gene_name'],) if candidate.get('gene_name') else (),
                      transcript_ids=(candidate['transcript_id'],), orf_ids=(candidate['orf_id'],),
                      orf_kind=model.get('orf_type'))
    return graph.finish()


def _fusion(value, reference):
    graph = _Snapshot(value, reference, value.get('provenance', {}).get('source'))
    observations = _rows(value['evidence']['observations'], 'Fusion observations')
    if value['provenance'].get('sample_id') != graph.sample or any(o['sample_id'] != graph.sample for o in observations):
        raise ValueError("Fusion RNA evidence scope does not match sample_id")
    support = value['evidence']['support']
    fragments = {(o['sample_id'], o['library_id'], o['fragment_id']) for o in observations}
    reads = {(o['sample_id'], o['library_id'], o['read_id']) for o in observations}
    if (type(support['reads']) is not int or type(support['fragments']) is not int
            or len(reads) != len(observations) or support['reads'] != len(reads)
            or support['fragments'] != len(fragments)):
        raise ValueError("Fusion RNA support disagrees with its original observations")
    event = graph.record('event', value['event_id'], {k: v for k, v in value.items() if k not in ('paths', 'evidence')},
                         status=value['status'])
    evidence = graph.record('evidence', 'original_rna', value['evidence'], support='rna_sequence',
                            scope={'sample': graph.sample, 'read_set': graph.source})
    graph.link(evidence, event, 'supports_rna_sequence')
    _path_translations(graph, event, value.get('paths'), origin='fusion')
    return graph.finish()


def _path_translations(graph, event, values, *, origin):
    """Map all annotated-frame products without choosing a transcript winner."""
    paths = _rows(values, 'RNA paths')
    _indexed(paths, 'path_id')
    for i, path in enumerate(paths):
        record = graph.record('path', path['path_id'], {k: v for k, v in path.items() if k != 'translations'},
                              status=path['frame_status'])
        graph.link(event, record, 'has_path', index=i)
        for j, row in enumerate(_rows(path.get('translations', []), 'RNA translations')):
            dna = _coding_sequence(path['sequence'], (row['translation_start'], row['translation_end']))
            if type(row['complete_5prime']) is not bool or row.get('translation_observed') is not False:
                raise ValueError("RNA candidates must retain unobserved translation and start status")
            if row['complete_5prime'] and dna[:3] not in standard_genetic_code.start_codons:
                raise ValueError("RNA start codon is unavailable")
            _translation(dna, row['amino_acids'], row['ends_with_stop_codon'], start=row['complete_5prime'])
            tid = path['path_id'] + ':translation:' + str(j)
            translation = graph.record('translation', tid, row, status=path['frame_status'])
            graph.link(record, translation, 'has_translation', index=j)
            transcripts = sorted(set(row['transcript_ids']) | {r['transcript_id'] for r in row.get('acceptor_frames', [])})
            graph.protein(translation, tid, row['amino_acids'], origin, start=row['complete_5prime'],
                          stop=row['ends_with_stop_codon'], transcript_ids=tuple(transcripts))


def _native_sv(value, reference):
    graph = _Snapshot(value, reference, value.get('source'))
    for row in value.get('cell_umi_evidence', {}).get('reads', []):
        scope = row.get('scope')
        if scope is None:
            continue
        if not isinstance(scope, Mapping):
            raise ValueError("Invalid SV RNA cell/UMI evidence scope")
        if any(scope.get(key) is not None and scope[key] != expected for key, expected in
               (('sample_id', graph.sample), ('source', graph.source))):
            raise ValueError("SV RNA cell/UMI evidence scope contradicts sample_id or source")
    evidence_keys = {'observations', 'original_records', 'record_evidence', 'cell_umi_evidence', 'read_lineage'}
    event = graph.record('event', value['event_id'],
                         {k: v for k, v in value.items() if k not in evidence_keys | {'paths'}},
                         status=value.get('status', 'status_unavailable'))
    evidence = graph.record('evidence', 'original_rna',
                            {k: value[k] for k in evidence_keys if k in value}, support='rna_sequence',
                            scope={'sample': graph.sample, 'read_set': graph.source})
    graph.link(evidence, event, 'supports_rna_sequence')
    _path_translations(graph, event, value.get('paths'), origin='sv_translation')
    exploratory = _sv(export_sv_rna_orfs(value), reference)
    graph.records.extend(exploratory.records)
    graph.links.extend(exploratory.links)
    graph.proteins.extend(exploratory)
    exported_event, = [record for record in exploratory.records if record.kind == 'event']
    graph.link(event, exported_event, 'has_exploratory_orf_export')
    return graph.finish()


def from_rna_result(result, *, reference):
    """Adapt native or offline RNA results without inferring protein expression.

    Parameters
    ----------
    result : dict
        An Isovar protein-hypotheses v2, SV ORF v4, SV reconstruction v5,
        fusion RNA v3 or variant RNA ORF v1 result. Native SV results retain
        annotated-frame products and use the same exploratory-ORF mapping as
        compact v4 exports. A JSON round trip uses exactly the same mapping.
    reference : orfdata.ReferenceContext
        Explicit genomic reference. Recorded builds and complete ORF model
        contexts must agree. Unrecorded species/annotation remain caller claims.

    Returns
    -------
    orfdata.Proteome
        Every retained protein hypothesis, plus source-qualified producer
        records, RNA evidence, placements and unresolved outcomes. No grouping
        by amino-acid sequence and no observed translation evidence are added.
    """
    if not isinstance(result, Mapping):
        raise ValueError("Expected an Isovar RNA result JSON object")
    if not isinstance(reference, ReferenceContext) or not reference.is_genomic:
        raise ValueError("A genomic ORFdata ReferenceContext is required")
    schema = result.get('schema')
    if schema == 'isovar.protein_hypotheses.v2':
        return from_protein_hypotheses(result, reference=reference)
    handlers = {'isovar.sv_rna_orfs.v4': _sv, 'isovar.variant_rna_orfs.v1': _variant,
                'isovar.fusion_rna.v3': _fusion, 'isovar.sv_rna_candidates.v5': _native_sv}
    if schema not in handlers:
        raise ValueError(f"Unsupported Isovar RNA result schema: {schema!r}")
    try:
        build = result['variant']['reference_name'] if schema == 'isovar.variant_rna_orfs.v1' else result['reference_name']
        if build != reference.genome_build:
            raise ValueError(f"Isovar build {build!r} conflicts with {reference.genome_build!r}")
        return handlers[schema](result, reference)
    except (KeyError, TypeError) as error:
        raise ValueError(f"Malformed {schema} result: {error}") from error
