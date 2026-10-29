"""One immutable library index must preserve complete original-template evidence."""

import json
from pathlib import Path

import pysam
import pytest

from isovar.read_collector import ReadCollector
from isovar.rna_library import prepare_rna_library
from .test_sv_rna import record, write_bam


def library_records():
    first = record('shared', '1', 100, '5S20M', 'T' * 5 + 'A' * 20, flag=65)
    mate = record('shared', '2', 15000, '25M', 'C' * 25, flag=129)
    alternative = record('shared', '2', 17000, '25M', 'G' * 25, flag=65 | 256)
    other_group = record('shared', '1', 100, '25M', 'T' * 25, group='b')
    unmapped = record('shared', '2', 0, None, 'ACGT' * 6, flag=129 | 4)
    background = record('background', '1', 5000, '25M', 'A' * 25)
    return [first, mate, alternative, other_group, unmapped, background, first]


def test_regional_query_keeps_distant_placements_mates_and_original_accounting(tmp_path):
    reads = library_records()
    path = write_bam(tmp_path / 'original.bam', reads)
    prepared = prepare_rna_library(path, tmp_path / 'indexes', read_collector=ReadCollector(min_mapping_quality=0))
    assert not list(prepared.directory.glob('*.bam'))
    assert prepared.statistics['records_scanned'] == 7
    assert prepared.statistics['unique_records'] == 6
    assert prepared.statistics['fragments_scanned'] == 3
    assert prepared.statistics['excluded_records'] == {'unmapped': 1}
    with prepared.open() as index:
        templates = list(index.templates([('1', 110, 111)]))
        assert len(templates) == 2
        observed = {read.to_string() for template in templates for read in template}
        assert observed == {read.to_string() for read in reads if read.query_name == 'shared' and not read.is_unmapped}
        full = index.template(('a', 'shared'), eligible_only=False)
        assert any(read.is_unmapped for records in full.values() for read in records)
        assert sum(len(records) for records in full.values()) == 4
        assert len(list(index.original_records())) == 7
        # Soft-clip sequence never expands the genomic interval into coverage.
        assert not list(index.templates([('1', 95, 100)]))
        assert ('a', 'shared', 64) in index
        assert ('b', 'shared', 64) not in index


def test_index_and_exhaustive_acquisition_agree_on_exact_cigar_blocks(tmp_path):
    reads = library_records() + [record('spliced', '1', 300, '10M1000N10M', 'A' * 20)]
    path = write_bam(tmp_path / 'original.bam', reads)
    regions = [('1', 310, 1310), ('2', 16999, 17001)]
    collector = ReadCollector(min_mapping_quality=0)
    expected = set()
    with pysam.AlignmentFile(path) as bam:
        for read in bam:
            if collector.alignment_filter_reason(read) is None and any(
                    read.reference_name == contig and start < hi and lo < end
                    for contig, start, end in regions for lo, hi in read.get_blocks()):
                expected.add((read.get_tag('RG'), read.query_name))
    with prepare_rna_library(path, tmp_path / 'indexes', read_collector=collector).open() as index:
        selected = list(index.templates(regions))
        exact = {(template[0].get_tag('RG'), template[0].query_name) for template in selected if any(
            read.reference_name == contig and start < hi and lo < end
            for read in template for contig, start, end in regions for lo, hi in read.get_blocks())}
        assert exact == expected
        assert any(template[0].query_name == 'spliced' for template in selected)  # conservative overselection


def test_cache_reuse_is_read_only_and_separates_eligibility_policy(tmp_path, monkeypatch):
    from isovar import rna_library
    path = write_bam(tmp_path / 'original.bam', library_records())
    cache = tmp_path / 'indexes'
    prepared = prepare_rna_library(path, cache)
    before = {p.name: p.read_bytes() for p in prepared.directory.iterdir()}
    original_build = rna_library._index_records

    def no_rebuild(*args):
        raise AssertionError('The unchanged library index must be reused')

    monkeypatch.setattr(rna_library, '_index_records', no_rebuild)
    reused = prepare_rna_library(path, cache)
    with prepared.open() as first, reused.open() as second:
        assert first[('a', 'shared', 64)][0].to_string() == second[('a', 'shared', 64)][0].to_string()
        assert list(first.templates([('1', 100, 110)]))
        with pytest.raises(Exception, match='readonly'):
            first.connection.execute('DELETE FROM records')
    assert before == {p.name: p.read_bytes() for p in prepared.directory.iterdir()}
    monkeypatch.setattr(rna_library, '_index_records', original_build)
    changed = prepare_rna_library(path, cache, read_collector=ReadCollector(use_secondary_alignments=False))
    assert changed.directory != prepared.directory
    assert changed.statistics['excluded_records']['secondary'] == 1


def test_changed_bytes_and_interruption_cannot_reuse_a_valid_index(tmp_path, monkeypatch):
    from isovar import rna_library
    path = Path(write_bam(tmp_path / 'original.bam', library_records()))
    prepared = prepare_rna_library(path, tmp_path / 'indexes')
    with pytest.raises(ValueError, match='changed during use'), prepared.open() as index:
        path.touch()
        index.template(('a', 'shared'))
    original = rna_library._index_records

    def interrupted(*args):
        original(*args)
        raise OSError('injected interruption')

    monkeypatch.setattr(rna_library, '_index_records', interrupted)
    scratch = tmp_path / 'interrupted'
    with pytest.raises(OSError, match='injected interruption'):
        prepare_rna_library(path, scratch)
    assert not list(scratch.iterdir())


def test_corrupted_cache_is_an_error_not_an_empty_library(tmp_path):
    path = write_bam(tmp_path / 'original.bam', library_records())
    prepared = prepare_rna_library(path, tmp_path / 'indexes')
    with (prepared.directory / 'records.sqlite').open('ab') as stream:
        stream.write(b'changed')
    with pytest.raises(ValueError, match='changed during use'):
        prepared.open()
    with pytest.raises(ValueError, match='digest differs'):
        prepare_rna_library(path, tmp_path / 'indexes')


def test_cached_statistics_are_pinned_with_the_index_bytes(tmp_path):
    path = write_bam(tmp_path / 'original.bam', library_records())
    prepared = prepare_rna_library(path, tmp_path / 'indexes')
    receipt_path = prepared.directory / 'receipt.json'
    receipt = json.loads(receipt_path.read_text())
    receipt['statistics']['records_scanned'] = 0
    receipt_path.write_text(json.dumps(receipt))
    with pytest.raises(ValueError, match='database metadata'):
        prepare_rna_library(path, tmp_path / 'indexes')


def test_custom_filter_requires_explicit_reproducible_identity(tmp_path):
    path = write_bam(tmp_path / 'original.bam', library_records())
    collector = ReadCollector(read_filter=lambda read: read.query_name != 'background')
    with pytest.raises(ValueError, match='read_filter_identity'):
        prepare_rna_library(path, tmp_path / 'indexes', read_collector=collector)
    prepared = prepare_rna_library(path, tmp_path / 'indexes', read_collector=collector, read_filter_identity='exclude-background.v1')
    assert prepared.statistics['excluded_records']['read_filter'] == 1


def test_collector_policy_is_checked_and_public_snapshots_cannot_mutate_counts(tmp_path):
    path = write_bam(tmp_path / 'original.bam', library_records())
    collector = ReadCollector()
    prepared = prepare_rna_library(path, tmp_path / 'indexes', read_collector=collector)
    prepared.statistics['records_scanned'] = 0
    prepared.identity['eligibility']['min_mapping_quality'] = 999
    with prepared.open() as index:
        index.statistics['records_scanned'] = 0
        assert index.statistics['records_scanned'] == 7
        index.validate_eligibility(collector)
        # Assembly/read-end choices are separate from original-record eligibility.
        collector.use_soft_clipped_bases = not collector.use_soft_clipped_bases
        index.validate_eligibility(collector)
        collector.min_mapping_quality += 1
        with pytest.raises(ValueError, match='eligibility differs'):
            index.validate_eligibility(collector)


def test_overridden_eligibility_requires_identity_and_preserves_its_exclusions(tmp_path):
    class CustomCollector(ReadCollector):
        def alignment_filter_reason(self, read):
            return 'custom_background' if read.query_name == 'background' else super().alignment_filter_reason(read)

    path = write_bam(tmp_path / 'original.bam', library_records())
    collector = CustomCollector()
    with pytest.raises(ValueError, match='read_filter_identity'):
        prepare_rna_library(path, tmp_path / 'indexes', read_collector=collector)
    prepared = prepare_rna_library(path, tmp_path / 'indexes', read_collector=collector, read_filter_identity='custom.v1')
    assert prepared.statistics['excluded_records']['custom_background'] == 1
    with prepared.open() as index:
        index.validate_eligibility(collector, read_filter_identity='custom.v1')
        with pytest.raises(ValueError, match='eligibility differs'):
            index.validate_eligibility(collector, read_filter_identity='custom.v2')


def test_actual_source_bytes_and_implementation_changes_get_new_indexes(tmp_path, monkeypatch):
    from isovar import rna_library
    path = write_bam(tmp_path / 'original.bam', library_records())
    cache = tmp_path / 'indexes'
    first = prepare_rna_library(path, cache)
    write_bam(path, library_records() + [record('added', '1', 800, '20M', 'T' * 20)])
    second = prepare_rna_library(path, cache)
    assert second.directory != first.directory
    assert second.statistics['records_scanned'] == 8
    monkeypatch.setattr(rna_library, 'software_identity', lambda: {'implementation': 'changed'})
    third = prepare_rna_library(path, cache)
    assert third.directory != second.directory


def test_source_exact_viral_decoy_alt_names_and_ambiguity_use_varcode(tmp_path):
    names = ['0', '1', 'CHR1', 'chrM', 'chrUn_KI270442v1', 'KI270442.2', 'HPV16', 'decoy-custom']
    header = pysam.AlignmentHeader.from_dict({'SQ': [{'SN': name, 'LN': 1000} for name in names]})
    reads = []
    for tid, name in enumerate(names):
        read = pysam.AlignedSegment(header)
        read.query_name, read.reference_id, read.reference_start = name, tid, 100
        read.query_sequence, read.cigarstring, read.mapping_quality = 'A' * 20, '20M', 60
        read.query_qualities = [30] * 20
        reads.append(read)
    path = write_bam(tmp_path / 'vocabulary.bam', reads, header)
    with prepare_rna_library(path, tmp_path / 'indexes', reference_name='GRCh38').open() as index:
        for query, expected in [('0', '0'), ('MT', 'chrM'), ('KI270442.1', 'chrUn_KI270442v1'),
                                ('KI270442.2', 'KI270442.2'), ('HPV16', 'HPV16'), ('decoy-custom', 'decoy-custom')]:
            template, = list(index.templates([(query, 100, 101)]))
            assert template[0].reference_name == expected
        with pytest.raises(ValueError):
            list(index.templates([('chr1', 100, 101)]))
