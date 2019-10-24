from isovar.common import normalize_base0_range_indices
from nose.tools import eq_

def normalize_base0_range_indices_valid_range():
    eq_(normalize_base0_range_indices(3, 7, 10),
        (3, 7))

def normalize_base0_range_indices_negative_indices():
    eq_(normalize_base0_range_indices(-3, -2, 10),
        (7, 8))

def normalize_base0_range_indices_start_None():
    eq_(normalize_base0_range_indices(None, -2, 10),
        (0, 8))
def normalize_base0_range_indices_end_None():
    eq_(normalize_base0_range_indices(-3, None, 10),
        (7, 10))
