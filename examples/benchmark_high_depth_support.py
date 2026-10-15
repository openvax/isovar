"""Reproduce #463 using Topiary's original MT-ND5 fixture, without sampling.

Run in a checkout with Isovar's dependencies installed:
    python -m examples.benchmark_high_depth_support TOPIARY_CHECKOUT output.json

Run the same script with PYTHONPATH pointing to a baseline checkout to compare.
Reference indexing and output fingerprints are outside the timed pipeline.
Use separate processes for comparable peak RSS; --profile adds cProfile overhead.
"""

import argparse
import cProfile
from hashlib import sha256
from importlib.metadata import version
import json
import logging
from pathlib import Path
import platform
import resource
import time

import osteosarc
from pyensembl import Genome
import pysam
from varcode import Variant

import isovar
from isovar.reference_coding_sequence_key import ReferenceCodingSequenceKey
from isovar.translation_key import TranslationKey


def fingerprint(value):
    return sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def fields(value, names):
    return {name: getattr(value, name) for name in names}


class ObservedCreator(isovar.ProteinSequenceCreator):
    """Retain existing outputs for comparison after the timed pipeline."""

    def all_pairs_translations(self, variant_sequences, reference_contexts):
        self.candidates = variant_sequences
        self.translations = super().all_pairs_translations(variant_sequences, reference_contexts)
        return self.translations


def evidence(result, creator):
    read_hashes = {}

    def read_hash(read):
        if id(read) not in read_hashes:
            values = fields(read, read.__slots__)
            ids = values['compatible_transcript_ids']
            values['compatible_transcript_ids'] = None if ids is None else sorted(ids)
            read_hashes[id(read)] = fingerprint(values)
        return read_hashes[id(read)]

    def reads_hash(reads):
        return fingerprint(sorted(read_hash(read) for read in reads))

    candidates = []
    candidate_hashes = {}
    for candidate in creator.candidates:
        ids = candidate.compatible_transcript_ids
        record = dict(prefix=candidate.prefix, alt=candidate.alt, suffix=candidate.suffix,
                      transcript_ids=None if ids is None else sorted(ids),
                      reads_sha256=reads_hash(candidate.reads),
                      placements_sha256=fingerprint(sorted(candidate._source_alignments.items())),
                      coverage=candidate.coverage().tolist())
        candidates.append(record)
        candidate_hashes[id(candidate)] = fingerprint(record)
    translations = []
    for translation in creator.translations:
        context = translation.reference_context
        translations.append(dict(
            **fields(translation, TranslationKey._fields),
            candidate=candidate_hashes[id(translation.untrimmed_variant_sequence)],
            orf=fields(translation.variant_orf, translation.variant_orf._fields),
            context=fields(context, ReferenceCodingSequenceKey._fields),
            transcripts=sorted(t.id for t in context.transcripts)))
    proteins = [dict(**fields(p, TranslationKey._fields), reads_sha256=reads_hash(p.supporting_reads),
                     reads=p.num_supporting_reads, fragments=p.num_supporting_fragments,
                     transcripts=sorted(p.transcript_ids)) for p in result.sorted_protein_sequences]
    return dict(candidates=candidates, translations=translations, proteins=proteins,
                counts={f'num_{category}_{unit}': getattr(result, f'num_{category}_{unit}')
                        for category in ('total', 'ref', 'alt', 'other') for unit in ('reads', 'fragments')},
                read_evidence={category: reads_hash(getattr(result.read_evidence, category + '_reads'))
                               for category in ('ref', 'alt', 'other')}, filters=result.filter_values)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('topiary', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--profile', type=Path)
    parser.add_argument('--variant-id', default='MT_ND5-chrM-12994')
    args = parser.parse_args()
    logging.disable(logging.CRITICAL)
    fixture = args.topiary / 'tests/data/osteosarc_all_variants'
    reference = fixture / 'reference'
    manifest = json.loads((reference / 'manifest.json').read_text())
    for name, receipt in manifest['files'].items():
        if osteosarc.digest(reference / name) != receipt['sha256']:
            raise ValueError('Reference checksum mismatch: ' + name)
    member = 'topiary/osteosarc_all_variants/source/t2-all-variant-regions.bam'
    bam_path = osteosarc.bundle_file('openvax-v1', member, offline=True)
    genome = Genome(reference_name='GRCh38-osteosarc-all-' + osteosarc.digest(reference / 'manifest.json')[:16],
                    annotation_name='osteosarc-ensembl-subset', annotation_version=87,
                    gtf_path_or_url=str(reference / 'reference.gtf.gz'),
                    transcript_fasta_paths_or_urls=[str(reference / 'reference.cdna.fa.gz')],
                    protein_fasta_paths_or_urls=[str(reference / 'reference.pep.fa.gz')],
                    copy_local_files_to_cache=True,
                    cache_directory_path=str(args.output.parent / 'reference-index'))
    genome.index()
    creator = ObservedCreator(protein_context_peptide_length=25, variant_sequence_assembly=True)
    creator.candidates, creator.translations = [], []
    collector = isovar.ReadCollector(use_secondary_alignments=False, use_duplicate_reads=False,
                                    min_mapping_quality=20, use_soft_clipped_bases=False,
                                    merge_overlapping_fragments=True)
    entry, = [r for r in json.loads((fixture / 'checked-inventory.json').read_text())['variants']
              if r['variant_id'] == args.variant_id and r['input_status'] == 'ready']
    contig = entry['chrom'].removeprefix('chr')
    variant = Variant('MT' if contig == 'M' else contig, entry['pos'], entry['ref'], entry['alt'], ensembl=genome)
    profile = cProfile.Profile() if args.profile else None
    with pysam.AlignmentFile(bam_path) as bam:
        wall, cpu = time.perf_counter(), time.process_time()
        if profile:
            profile.enable()
        result, = isovar.run_isovar([variant], bam, read_collector=collector, protein_sequence_creator=creator)
        if profile:
            profile.disable()
        seconds, cpu_seconds = time.perf_counter() - wall, time.process_time() - cpu
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    output = evidence(result, creator)
    source = Path(isovar.__file__).parent
    report = dict(isovar_version=isovar.__version__, isovar_path=isovar.__file__, python=platform.python_version(),
                  code_sha256={name: osteosarc.digest(source / name) for name in
                               ('read_identity.py', 'assembly.py', 'variant_sequence.py', 'variant_sequence_helpers.py')},
                  variant=entry,
                  dependencies={name: version(name) for name in ('varcode', 'pyensembl', 'osteosarc', 'numpy', 'pysam')},
                  benchmark_sha256=osteosarc.digest(Path(__file__)), bundle='openvax-v1', member=member,
                  bam_sha256=osteosarc.digest(Path(bam_path)), reference_manifest=manifest,
                  seconds=seconds, cpu_seconds=cpu_seconds, profiled=profile is not None,
                  peak_process_rss_bytes=rss if platform.system() == 'Darwin' else rss * 1024,
                  evidence_sha256=fingerprint(output), evidence=output)
    with args.output.open('x') as handle:
        json.dump(report, handle, indent=2)
        handle.write('\n')
    if profile:
        profile.dump_stats(args.profile)
    print(json.dumps({k: v for k, v in report.items() if k not in ('evidence', 'reference_manifest')}))


if __name__ == '__main__':
    main()
