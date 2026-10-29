"""Regional acquisition filters locally and never publishes an incomplete receipt."""

from dataclasses import asdict
from io import BytesIO

import pysam
import pytest
from varcode import Variant

from reports.matched_normal import acquire_sid_normal as acquisition
from isovar.sv_rna_audit import file_sha256
from .test_normal_sources import source_file


class Response(BytesIO):
    headers = {'ETag': '"public-object-identity"'}


def test_acquisition_hashes_full_stream_and_sends_no_genomic_regions(tmp_path, monkeypatch):
    original = source_file(tmp_path)
    data = open(original.path, 'rb').read()
    descriptor = asdict(original)
    descriptor['provenance'] = dict(url='https://example.invalid/public.vcf.gz',
        source_size=len(data), http_identity=Response.headers)
    requests = []

    def open_source(request, timeout):
        requests.append(request)
        return Response(data)

    monkeypatch.setattr(acquisition, 'urlopen', open_source)
    output = tmp_path / 'acquired'
    output.mkdir()
    variants = dict(sha256='f' * 64, values=[Variant('1', 105, 'A', 'G')])
    result = acquisition.acquire(descriptor, variants, 1, output)
    request, = requests
    assert request.full_url == descriptor['provenance']['url']
    assert request.data is None
    assert request.header_items() == [('If-match', Response.headers['ETag'])]
    assert result['source']['provenance']['original_file']['sha256'] == file_sha256(original.path)
    with pysam.VariantFile(result['source']['path']) as reader:
        record, = list(reader)
        assert record.samples['EXPLICIT.sample-column']['GT'] == (0, 2)
        assert record.alts == ('G', 'T')  # Acquisition does no genotype filtering.
    assert result['source']['acquired_regions'] == (('chr1', 103, 106),)
    acquisition.acquire(descriptor, variants, 1, output)
    assert len(requests) == 1
    receipt_path = next(output.glob('*.acquisition.json.gz'))
    digest = file_sha256(receipt_path)
    new_policy = dict(descriptor, policy=dict(descriptor['policy'], minimum_depth=8))
    reused = acquisition.acquire(new_policy, variants, 1, output)
    assert reused['source']['policy']['minimum_depth'] == 8
    assert len(requests) == 1 and file_sha256(receipt_path) == digest
    with pytest.raises(ValueError, match='scope changed'):
        acquisition.acquire(descriptor, variants, 2, output)


def test_incomplete_stream_never_publishes_an_acquisition_receipt(tmp_path, monkeypatch):
    original = source_file(tmp_path)
    data = open(original.path, 'rb').read()
    descriptor = asdict(original)
    descriptor['provenance'] = dict(url='https://example.invalid/public.vcf.gz',
        source_size=len(data) + 1, http_identity=Response.headers)
    monkeypatch.setattr(acquisition, 'urlopen', lambda *args, **kwargs: Response(data))
    output = tmp_path / 'incomplete'
    output.mkdir()
    with pytest.raises(ValueError, match='incomplete'):
        acquisition.acquire(descriptor, dict(sha256='f' * 64, values=[Variant('1', 105, 'A', 'G')]), 1, output, attempts=1)
    assert not list(output.glob('*.acquisition.json.gz'))
    assert not list(output.glob('*.vcf.gz'))


def test_acquisition_includes_distant_background_regions_and_invalidates_changed_scope(tmp_path, monkeypatch):
    original = source_file(tmp_path)
    data = open(original.path, 'rb').read()
    descriptor = asdict(original)
    descriptor['provenance'] = dict(url='https://example.invalid/public.vcf.gz',
        source_size=len(data), http_identity=Response.headers)
    monkeypatch.setattr(acquisition, 'urlopen', lambda *args, **kwargs: Response(data))
    output = tmp_path / 'distant'
    output.mkdir()
    variants = dict(sha256='f' * 64, values=[Variant('1', 800, 'A', 'G')])
    result = acquisition.acquire(descriptor, variants, 0, output, additional_regions=[('1', 100, 110)])
    assert result['counts']['retained_records'] == 1
    assert result['source']['acquired_regions'] == (('chr1', 100, 110), ('chr1', 799, 800))
    with pytest.raises(ValueError, match='scope changed'):
        acquisition.acquire(descriptor, variants, 0, output, additional_regions=[('1', 200, 210)])


def test_transient_transport_timeout_retries_the_same_public_identity(tmp_path, monkeypatch):
    original = source_file(tmp_path)
    data = open(original.path, 'rb').read()
    descriptor = asdict(original)
    descriptor['provenance'] = dict(url='https://example.invalid/public.vcf.gz',
        source_size=len(data), http_identity=Response.headers)
    requests, delays = [], []

    def open_source(request, timeout):
        requests.append((request.full_url, request.header_items(), request.data))
        if len(requests) == 1:
            raise TimeoutError('injected public transport timeout')
        return Response(data)

    monkeypatch.setattr(acquisition, 'urlopen', open_source)
    monkeypatch.setattr(acquisition, 'sleep', delays.append)
    output = tmp_path / 'retry'
    output.mkdir()
    result = acquisition.acquire(descriptor, dict(sha256='f' * 64, values=[Variant('1', 105, 'A', 'G')]), 1, output)
    assert requests[0] == requests[1]
    assert delays == [2]
    assert result['source']['provenance']['original_file']['sha256'] == file_sha256(original.path)
