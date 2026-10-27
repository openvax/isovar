"""Sequence display preserves codon anchors and shows deleted residues explicitly."""

from isovar.cohort_pdf import aligned_columns


def test_deletion_in_repeat_keeps_declared_anchor_and_neutral_flanks():
    ref = 'LAERKYSDEKKRVDELQQQVKLYKLNLESSKQE'
    rna = 'LAERKYSDEKRVDELQQQVKLYKLNLESSKQE'
    columns = aligned_columns(ref, rna, [9, 9], allele_length_delta=3)
    deleted = [i for i, c in enumerate(columns) if c['focal_deletion']]
    assert deleted == [9]
    assert columns[9]['reference'] == 'K'
    assert columns[9]['protein'] == '-'
    assert all(c['reference'] == c['protein'] for c in columns if not c['focal_deletion'])
    assert ''.join(c['protein'] for c in columns if c['position'] is not None) == rna


def test_frameshift_is_not_aligned_back_to_a_spurious_matching_tail():
    columns = aligned_columns('MKTLLAG', 'MKAGLLA', [2, 7], allele_length_delta=1, frameshift=True)
    assert [c['position'] for c in columns] == list(range(7))
    assert [c['reference_position'] for c in columns] == list(range(7))
    assert not any(c['focal_deletion'] for c in columns)


def test_truncated_reference_does_not_invent_matching_stop():
    columns = aligned_columns('MKAAALLL', 'MK', [2, 2], stop=True, reference_stop=True)
    assert columns[-1]['protein'] == '*'
    assert columns[-1]['reference'] != '*'
