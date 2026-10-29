"""The replay preparation changes declared coordinates, never RNA payloads."""

from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pysam
import pytest


def module():
    spec = spec_from_file_location('prepare_ddx17', Path(__file__).parents[1] / 'examples/prepare_ddx17_rna_paths.py')
    loaded = module_from_spec(spec)
    spec.loader.exec_module(loaded)
    return loaded


def test_projection_keeps_original_sam_fields_and_local_scope(tmp_path):
    prepare = module()
    audit, output = tmp_path / 'audit', tmp_path / 'output'
    audit.mkdir()
    output.mkdir()
    (audit / 'reference.fa').write_text('>local\n' + 'A' * 50000 + '\n')
    header = pysam.AlignmentHeader.from_dict({'SQ': [{'SN': prepare.LOCAL_CONTIG, 'LN': 50000}]})
    original = pysam.AlignedSegment(header)
    original.query_name = 'one-original'
    original.query_sequence = 'ACGTA'
    original.query_qualities = [20, 21, 22, 23, 24]
    original.reference_id, original.reference_start = 0, 100
    original.cigarstring, original.mapping_quality = '5M', 255
    original.set_tag('NH', 1)
    original.flag = 99
    original.next_reference_id, original.next_reference_start = 0, 200
    with pysam.AlignmentFile(audit / 'expanded.star.bam', 'wb', header=header) as bam:
        bam.write(original)
    unchanged = (audit / 'expanded.star.bam').read_bytes()

    class Reference:
        def fetch(self, contig, start, end):
            assert contig == 'chr22' and (start, end) == (38470000, 38520000)
            return 'A' * 50000

        def get_reference_length(self, contig):
            return 50818468

    path, count = prepare.project(audit, Reference(), 'chr22', output)
    assert count == 1
    assert (audit / 'expanded.star.bam').read_bytes() == unchanged
    with pysam.AlignmentFile(path) as bam:
        projected = next(bam)
        assert projected.reference_start == original.reference_start + prepare.OFFSET
        assert projected.next_reference_start == original.next_reference_start + prepare.OFFSET
        assert projected.query_sequence == original.query_sequence
        assert projected.query_qualities == original.query_qualities
        assert projected.cigarstring == original.cigarstring and projected.tags == original.tags
        assert projected.flag == original.flag
        assert 'LOCAL' in bam.header.to_dict()['CO'][0]

    original.set_tag('SA', 'chr22,1,+,5M,60,0;')
    with pysam.AlignmentFile(audit / 'expanded.star.bam', 'wb', header=header) as bam:
        bam.write(original)
    rejected = tmp_path / 'rejected'
    rejected.mkdir()
    with pytest.raises(ValueError, match='Unexpected tag'):
        prepare.project(audit, Reference(), 'chr22', rejected)
