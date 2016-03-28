# Copyright (c) 2016. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Helper functions for normalizing and working with genomic variants
"""

from __future__ import print_function, division, absolute_import

from .logging import create_logger

logger = create_logger(__name__)

def trim_variant_fields(location, ref, alt):
    """
    Trims common prefixes from the ref and alt sequences

    Parameters
    ----------
    location : int
        Position (starting from 1) on some chromosome

    ref : str
        Reference nucleotides

    alt : str
        Alternate (mutant) nucleotide

    Returns adjusted triplet (location, ref, alt)
    """
    if len(alt) > 0 and ref.startswith(alt):
        # if alt is a prefix of the ref sequence then we actually have a
        # deletion like:
        #   g.10 GTT > GT
        # which can be trimmed to
        #   g.12 'T'>''
        ref = ref[len(alt):]
        location += len(alt)
        alt = ""
    if len(ref) > 0 and alt.startswith(ref):
        # if ref sequence is a prefix of the alt sequence then we actually have
        # an insertion like:
        #   g.10 GT>GTT
        # which can be trimmed to
        #   g.11 ''>'T'
        # Note that we are selecting the position *before* the insertion
        # (as an arbitrary convention)
        alt = alt[len(ref):]
        location += len(ref) - 1
        ref = ""
    return location, ref, alt


def trim_variant(variant):
    """
    Parameters
    ----------
    variant : varcode.Variant

    Returns trimmed triplet (location, ref, alt)
    """
    return trim_variant_fields(variant.start, variant.ref, variant.alt)


def base0_interval_for_variant_fields(base1_location, ref, alt):
    """
    Inteval of interbase offsets of the affected reference positions for a
    particular variant.

    Parameters
    ----------
    base1_location : int
        First reference nucleotide of variant or, for insertions, the base
        before the insertion.

    ref : str
        Reference nucleotides

    alt : str
        Alternative nucleotides
    """
    if len(ref) == 0:
        # if the variant is an insertion then we need to check to make sure
        # both sides of the insertion are matches
        base0_start = base1_location - 1
        base0_end = base1_location + 1
    else:
        # substitution or deletion
        base0_start = base1_location - 1
        base0_end = base0_start + len(ref)
    return base0_start, base0_end


def base0_interval_for_variant(variant):
    """
    Parameters
    ----------
    variant : varcode.Variant

    Returns triplet of (base1_location, ref, alt)
    """
    base1_location, ref, alt = trim_variant(variant)
    print(base1_location, ref, alt)
    return base0_interval_for_variant_fields(
        base1_location=base1_location,
        ref=ref,
        alt=alt)


def interbase_range_affected_by_variant_on_transcript(variant, transcript):
    """
    Convert from a variant's position in global genomic coordinates on the
    forward strand to an interval of interbase offsets on a particular
    transcript's mRNA.

    Parameters
    ----------
    variant : varcode.Variant

    transcript : pyensembl.Transcript

    Assumes that the transcript overlaps the variant.

    Returns (start, end) tuple of offsets into the transcript's cDNA sequence
    which indicates which bases in the reference sequence are affected by a
    variant.

    Example:
        The insertion of "TTT" into the middle of an exon would result in an
        offset pair such as (100,100) since no reference bases are changed
        or deleted by an insertion.

        On the other hand, deletion the preceding "CGG" at that same locus could
        result in an offset pair such as (97, 100)
    """
    if variant.is_insertion:
        if transcript.strand == "+":
            # base-1 position of an insertion is the genomic nucleotide
            # before any inserted mutant nucleotides
            start_offset = [transcript.spliced_offset(variant.start)]
        else:
            # assuming that this transcript was considered to overlap
            # with the variant since the insertion happens inside
            # one of its exons (rather than simply immediately before
            # of after)
            start_offset = [transcript.spliced_offset(variant.start + 1)]
        # an insertion happens *between* two reference bases
        # so the start:end offsets coincide
        end_offset = start_offset
    else:
        # reference bases affected by substitution or deletion defined by
        # range starting at first affected base
        offsets = []
        assert len(variant.ref) > 0
        for dna_pos in range(variant.start, variant.start + len(variant.ref)):
            try:
                offsets.append(transcript.spliced_offset(dna_pos))
            except ValueError:
                logger.info(
                    "Couldn't find position %d from %s on exons of %s" % (
                        dna_pos,
                        variant,
                        transcript))
        if len(offsets) == 0:
            raise ValueError(
                "Couldn't find any exonic reference bases affected by %s on %s" % (
                    variant,
                    transcript))
        start_offset = min(offsets)
        end_offset = max(offsets) + 1
    return (start_offset, end_offset)
