from .value_object import ValueObject
from .variant_helpers import base0_interval_for_variant

class Interval(ValueObject):
    """
    Representation of a half-open interval on a chromosome with helpers
    for converting to other coordinate systems.
    """
    __slots__ = [
        "chromosome",
        "base0_start_inclusive",
        "base0_end_exclusive",
    ]

    def __len__(self):
        return self.base0_end_exclusive - self.base0_start_inclusive

    @property
    def base0_half_open_coordinates(self):
        return (self.base0_start_inclusive,  self.base0_end_exclusive)

    @property
    def base1_half_open_coordinates(self):
        return (self.base1_start_inclusive, self.base1_end_exclusive)

    @property
    def base0_open_coorindates(self):
        if self.base0_start_inclusive == self.base0_end_exclusive:
            return (self.base0_start_inclusive, self.base0_start_inclusive + 1)
        return (self.base0_start_exclusive, self.base0_end_exclusive)

    @property
    def base1_open_coorindates(self):
        return (self.base1_start_exclusive, self.base1_end_exclusive)

    @property
    def base0_closed_coorindates(self):
        return (self.base0_start_exclusive, self.base0_end_exclusive)

    @property
    def base1_closed_coorindates(self):
        return (self.base1_start_exclusive, self.base1_end_exclusive)


    @property
    def base0_start_exclusive(self):
        return self.base0_position_before

    @property
    def base0_position_before(self):
        return self.base0_start_inclusive - 1

    @property
    def base0_position_after(self):
        return self.base0_end_exclusive

    @property
    def base1_position_before(self):
        return self.base0_start_inclusive

    @property
    def base1_position_after(self):
        return self.base0_end_exclusive + 1

    @property
    def base1_start_inclusive(self):
        return self.base0_start_inclusive + 1

    @property
    def base1_end_exclusive(self):
        return self.base0_end_exclusive + 1

    @classmethod
    def from_variant(cls, variant):
        base0_start_inclusive, base0_end_exclusive = base0_interval_for_variant(variant)
        return cls(
            chromosome=variant.contig,
            base0_start_inclusive=base0_start_inclusive,
            base0_end_exclusive=base0_end_exclusive)
