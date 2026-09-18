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
Gathered all the default function parameters in a single module, so that these
values can be easily shared between modules and also between commandline
arguments of different scripts.
"""

from collections import OrderedDict

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

# Merge overlapping mates before assembly, retaining raw alignment counts.
MERGE_OVERLAPPING_FRAGMENTS = True

# htslib threads used to decompress BAM/CRAM input
NUM_RNA_DECOMPRESSION_THREADS = 1

# Context on each side of a variant for the reference-context inspection tool.
# Protein creation derives its context size from the requested cDNA length.
REFERENCE_CONTEXT_SIZE = 30

# number of nucleotides to extract from RNAseq reads around each variant
VARIANT_SEQUENCE_LENGTH = 90

# minimum number of reads supporting each nucleotide of a
# variant coding sequence
MIN_VARIANT_SEQUENCE_COVERAGE = 2

# Supplied fusion transcripts are an additive, direct-junction evidence path.
FUSION_PEPTIDE_LENGTHS = (8, 9, 10, 11)
MIN_FUSION_FRAGMENTS = 2

# number of nucleotides shared between reference and variant sequence
# before variant for reference contexts used to establish ORF
MIN_TRANSCRIPT_PREFIX_LENGTH = 10

# maximum number of mismatching nucleotides between reference and variant
# prefix sequences
MAX_REFERENCE_TRANSCRIPT_MISMATCHES = 2

# whether to include mismatches after a variant locus toward the
# MAX_REFERENCE_TRANSCRIPT_MISMATCHES count
COUNT_MISMATCHES_AFTER_VARIANT = False

# A centered single-residue mutation with 24 residues on either side covers
# every mutation-containing 25mer (2 * 25 - 1 = 49 residues).
PROTEIN_CONTEXT_PEPTIDE_LENGTH = 25


def protein_sequence_length_for_peptide_length(peptide_length):
    """Context covering every peptide window around a centered changed residue."""
    return 2 * peptide_length - 1


PROTEIN_SEQUENCE_LENGTH = protein_sequence_length_for_peptide_length(PROTEIN_CONTEXT_PEPTIDE_LENGTH)
PROTEIN_SEQUENCE_PREFERENCE = "balanced"
MIN_PROTEIN_SEQUENCE_SUPPORT_FRACTION = 0.85

# number of protein sequences we want to return per variant
MAX_PROTEIN_SEQUENCES_PER_VARIANT = 1

# run overlap assembly algorithm to construct variant
# sequences from multiple reads which only partially
# overlap (rather than fully spanning a coding sequence)
VARIANT_SEQUENCE_ASSEMBLY = True

# Only merge variant cDNA sequences which at least share
# this number of nucleotides. Should be sufficiently high
# to minimize false assembly of isoforms which don't
# actually exist.
MIN_VARIANT_SEQUENCE_ASSEMBLY_OVERLAP_SIZE = 30

# include sequences from reads if they were clipped
# by the aligner
USE_SOFT_CLIPPED_BASES = False

# Missing QUAL is unknown confidence, not evidence that the sequence is bad.
USE_READS_WITHOUT_BASE_QUALITIES = True

# End annotation and derived trimming are opt-in; raw alignments never change.
INFER_READ_ENDS = False
TRIM_ADAPTERS = False
TRIM_POLY_A = False
READ_END_PROFILE = None
READ_END_WINDOW = 200
MIN_ADAPTER_OVERLAP = 12
MAX_ADAPTER_ERROR_RATE = 0.15
MIN_POLY_A_LENGTH = 12
MAX_POLY_A_ERROR_RATE = 0.1

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


DEFAULT_FILTER_THRESHOLDS = OrderedDict([
    ("min_num_alt_reads", MIN_NUM_RNA_ALT_READS),
    ("min_num_alt_fragments", MIN_NUM_RNA_ALT_FRAGMENTS),
    ("min_fraction_alt_reads", MIN_FRACTION_RNA_ALT_READS),
    ("min_fraction_alt_fragments", MIN_FRACTION_RNA_ALT_FRAGMENTS),
    ("max_num_ref_reads", MAX_NUM_RNA_REF_READS),
    ("max_num_ref_fragments", MAX_NUM_RNA_REF_FRAGMENTS),
    ("max_fraction_ref_reads", MAX_FRACTION_RNA_REF_READS),
    ("max_fraction_ref_fragments", MAX_FRACTION_RNA_REF_FRAGMENTS),
    ("max_num_other_reads", MAX_NUM_RNA_OTHER_READS),
    ("max_num_other_fragments", MAX_NUM_RNA_OTHER_FRAGMENTS),
    ("max_fraction_other_reads", MAX_FRACTION_RNA_OTHER_READS),
    ("max_fraction_other_fragments", MAX_FRACTION_RNA_OTHER_FRAGMENTS),
    ("min_ratio_alt_to_other_fragments", MIN_RATIO_RNA_ALT_TO_OTHER_FRAGMENTS),
])

DEFAULT_FILTER_FLAGS = [
    "predicted_effect_modifies_protein_sequence",
    "has_mutant_protein_sequence_from_rna",
    "protein_sequence_contains_mutation",
]

# Static mutation-evidence figures; shared by the plotting API and CLI.
PLOT_COMPARE_ASSEMBLY = False
PLOT_ALL_PROTEINS = False
PLOT_PROTEIN_ROWS_PER_PAGE = 10
PLOT_DPI = 600
PLOT_OVERVIEW_DPI = 300
PLOT_WIDTH = 16
PLOT_MAX_ROWS = 20
PLOT_VIEW = "all"
PLOT_VIEWS = ("all", "protein", "coverage", "reads", "assembly", "transcripts")
PLOT_OUTPUT_DIRECTORY = "isovar-figures"
