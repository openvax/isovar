"""Large acquisition receipts stay readable and may span PDF pages."""

from copy import deepcopy

from pypdf import PdfReader

from isovar.cohort_pdf import _provenance_rows, write_cohort_reports


def test_structured_provenance_preserves_fields_without_mutation():
    value = dict(source_url='https://example.org/source', interval=['chr22', 10, 20],
        files={'original.bam': {'path': '/original/rna.bam', 'sha256': 'a' * 64}},
        statuses=[{'state': 'unavailable'}, {'state': 'uncalled'}], missing=None, sequence='A' * 64)
    original = deepcopy(value)
    rows = dict(_provenance_rows(value))
    assert rows['source_url'] == value['source_url']
    assert rows['interval'] == 'chr22, 10, 20'
    assert rows['files / original.bam / path'] == '/original/rna.bam'
    assert rows['files / original.bam / sha256'] == 'a' * 16 + ' (SHA-256 prefix)'
    assert rows['statuses / 1 / state'] == 'unavailable'
    assert rows['statuses / 2 / state'] == 'uncalled'
    assert rows['missing'] == 'Not supplied'
    assert rows['sequence'] == 'A' * 64
    assert value == original


def test_large_provenance_table_splits_without_losing_sources(tmp_path):
    provenance = dict(source_url='https://example.org/original-RNA', scope='Regional RNA; no absence conclusion',
        long_description='Declared original acquisition. ' * 500 + 'Last retained sentence.',
        files={f'original-{index:03d}.bam': dict(path=f'/inputs/original-{index:03d}.bam', sha256=f'{index:064x}')
               for index in range(60)})
    source = dict(source_id='original', sample_id='sample', library_id='library', scope='Regional RNA',
                  provenance=provenance)
    manifest = dict(samples=[dict(sample_id='sample')], sources=[source], files={})
    summary = dict(sources=[source], small_variants=[], additional=[])
    path, = write_cohort_reports(manifest, summary, tmp_path)
    reader = PdfReader(path)
    assert len(reader.pages) > 2
    text = '\n'.join(page.extract_text() for page in reader.pages)
    assert 'https://example.org/original-RNA' in text
    assert 'Regional RNA; no absence conclusion' in text
    for index in range(60):
        assert f'/inputs/original-{index:03d}.bam' in text
    assert 'SHA-256 prefix' in text
    assert 'Last retained sentence.' in text
    assert '"source_url"' not in text
