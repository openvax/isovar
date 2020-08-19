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

from __future__ import print_function, division, absolute_import

"""
Gathered all the default function parameters in a single module, so that these
values can be easily shared between modules and also between commandline
arguments of different scripts.
"""

# lowest mapping quality (MAPQ) value to allow for RNAseq reads
# Rationale for a default value of 1:
#   RNA aligners such as STAR typically reserve MAPQ=255 for unique alignments
#   and then have some heuristic for decreasing MAPQ as a function of the number
#   of locations to which a read can map.
#   Continuing to use STAR as an example, here's a guide to its MAPQ values:
#       255 = uniquely mapped reads
#       3 = read maps to 2 locations
#       2 = read maps to 3 locations
#       1 = reads maps to 4-9 locations
#       0 = reads maps to 10 or more locations
#   (from http://seqanswers.com/forums/archive/index.php/t-27470.html)
#   Since there's no bound on just how bad MAPQ=0 really can be, it seems
#   wises to exclude those reads but still allow 9 candidate locations
#   per read.
#
#   The exact numbers differ with other aligners but the basic principle is the
#   same: MAPQ=0 might have an extremely large number of mapping locations
#   (probably due to low complexity of the sequence) but anything with MAPQ=1
#   is probably useful.
MIN_READ_MAPPING_QUALITY = 1

# use a read even if it's been marked as a duplicate?
USE_DUPLICATE_READS = False

# use a read even at a location that isn't its primary alignment?
USE_SECONDARY_ALIGNMENTS = True

# number of nucleotides to extract from RNAseq reads around each variant
VARIANT_SEQUENCE_LENGTH = 90

# minimum number of reads supporting each nucleotide of a
# variant coding sequence
MIN_VARIANT_SEQUENCE_COVERAGE = 2

# number of nucleotides shared between reference and variant sequence
# before variant for reference contexts used to establish ORF
MIN_TRANSCRIPT_PREFIX_LENGTH = 10

# maximum number of mismatching nucleotides between reference and variant
# prefix sequences
MAX_REFERENCE_TRANSCRIPT_MISMATCHES = 2

# whether to include mismatches after a variant locus toward the
# MAX_REFERENCE_TRANSCRIPT_MISMATCHES count
COUNT_MISMATCHES_AFTER_VARIANT = False

# number of amino acids / codons we're trying to translate
PROTEIN_SEQUENCE_LENGTH = 20

# number of protein sequences we want to return per variant
MAX_PROTEIN_SEQUENCES_PER_VARIANT = 1

# run overlap assembly algorithm to construct variant
# sequences from multiple reads which only partially
# overlap (rather than fully spanning a coding sequence)
VARIANT_SEQUENCE_ASSEMBLY = False

# Only merge variant cDNA sequences which at least share
# this number of nucleotides. Should be sufficiently high
# to minimize false assembly of isoforms which don't
# actually exist.
MIN_VARIANT_SEQUENCE_ASSEMBLY_OVERLAP_SIZE = 30

# include sequences from reads if they were clipped
# by the aligner
USE_SOFT_CLIPPED_BASES = False

# minimum number of RNA reads supporting a variant allele
MIN_NUM_RNA_ALT_READS = 3

# minimum number of total RNA fragments supporting a variant allele,
# differs from MIN_NUM_ALT_READS for paired end sequencing but is the same
# for single-end sequencing
MIN_NUM_RNA_ALT_FRAGMENTS = 2

# minimum ratio of # alt reads / # total overlapping reads
MIN_FRACTION_RNA_ALT_READS = 0.005  # (at least e.g. 3 in 600)

# minimum ratio of # alt fragments / # total overlapping fragments
MIN_FRACTION_RNA_ALT_FRAGMENTS = 0.005  # (at least e.g. 3 in 600)

# maximum number of RNA reads supporting a reference allele
MAX_NUM_RNA_REF_READS = 10 ** 9

# maximum number of total RNA fragments supporting a reference allele
MAX_NUM_RNA_REF_FRAGMENTS = 10 ** 9

# minimum ratio of # ref reads / # total overlapping reads
MAX_FRACTION_RNA_REF_READS = 1.0

# minimum ratio of # ref fragments / # total overlapping fragments
MAX_FRACTION_RNA_REF_FRAGMENTS = 1.0

# maximum number of RNA reads supporting a reference allele
MAX_NUM_RNA_OTHER_READS = 10 ** 9

# maximum number of total RNA fragments supporting a reference allele
MAX_NUM_RNA_OTHER_FRAGMENTS = 10 ** 9

# minimum ratio of # other (non-ref/non-alt) reads / # total overlapping reads
MAX_FRACTION_RNA_OTHER_READS = 0.5

# minimum ratio of # other fragments (non-ref/non-alt) / # total overlapping fragments
MAX_FRACTION_RNA_OTHER_FRAGMENTS = 0.5

# At loci where there is RNA support for both the alt allele and other
# non-reference alleles, we want the number of reads supporting the alt
# to be at least this many times greater than the total counts for the
# third and fourth alleles.
MIN_RATIO_RNA_ALT_TO_OTHER_FRAGMENTS = 3.0

# number of RNA fragments shared between two assembled protein sequences
# before we say that their variants are phased
MIN_SHARED_FRAGMENTS_FOR_PHASING = 2