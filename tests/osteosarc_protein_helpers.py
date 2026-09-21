"""Independent protein oracle: raw GTF coordinates + FASTA + NCBI table 1.

No Isovar or Varcode code is used to map, edit or translate the reference.
This deliberately supports only the six pinned single-exon SNVs/deletions,
not arbitrary variant normalization or transcript consequence annotation.
"""

from collections import defaultdict
import gzip
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import re

from .data.osteosarc.protein_references import fasta_records
from .reference_identity import reference_dataset_identity


REFERENCE = Path(__file__).parent / "data" / "osteosarc" / "protein_reference"
# https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1
CODE = dict(zip(("".join(b) for b in product("TCAG", repeat=3)),
                "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))


def translate(sequence):
    amino_acids = []
    for i in range(0, len(sequence) - 2, 3):
        amino_acid = CODE.get(sequence[i:i + 3], "X")
        if amino_acid == "*":
            return "".join(amino_acids), True
        amino_acids.append(amino_acid)
    return "".join(amino_acids), False


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def transcript_offset(exons, strand, position):
    offset = 0
    for start, end in sorted(exons, reverse=strand == "-"):
        if start <= position <= end:
            return offset + (position - start if strand == "+" else end - position)
        offset += end - start + 1
    raise ValueError(f"Position {position} is outside the pinned transcript exons")


def load_references(directory=REFERENCE):
    manifest = json.loads((directory / "protein_reference_manifest.json").read_text())
    for name, metadata in manifest["files"].items():
        data = (directory / name).read_bytes()
        if (sha256(data).hexdigest() != metadata["subset_sha256"]
                or sha256(gzip.decompress(data)).hexdigest() != metadata["uncompressed_sha256"]):
            raise ValueError(f"Protein reference checksum mismatch: {name}")
    cdnas = {h[1:].split(".")[0]: (h, s) for h, s in fasta_records(directory / "reference.cdna.fa.gz")}
    proteins = {re.search(r"transcript:(ENST\d+)", h)[1]: (h, s)
                for h, s in fasta_records(directory / "reference.pep.fa.gz")}
    features = defaultdict(lambda: defaultdict(list))
    with gzip.open(directory / "reference.gtf.gz", "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            transcript = re.search(r'transcript_id "([^"]+)"', fields[8])
            if transcript:
                features[transcript[1]][fields[2]].append(fields)
    result = {}
    for gene, tid in manifest["transcripts"].items():
        rows = features[tid]
        strand = rows["transcript"][0][6]
        exons = [(int(r[3]), int(r[4])) for r in rows["exon"]]
        start_bases = [p for r in rows["start_codon"] for p in range(int(r[3]), int(r[4]) + 1)]
        cds_start = min(transcript_offset(exons, strand, p) for p in start_bases)
        header, cdna = cdnas[tid]
        protein_header, protein = proteins[tid]
        assert sum(e - s + 1 for s, e in exons) == len(cdna)
        assert translate(cdna[cds_start:]) == (protein, True), tid
        result[gene] = {
            "transcript_id": tid, "transcript_version": header[1:].split()[0],
            "protein_version": protein_header[1:].split()[0],
            "contig": rows["transcript"][0][0], "strand": strand, "exons": exons,
            "cds_start": cds_start, "cdna": cdna, "protein": protein,
        }
    return result


def expected_variant(record, reference):
    start, ref, alt = int(record["pos"]), record["ref"], record["alt"]
    assert record["chrom"].removeprefix("chr") == reference["contig"]
    if len(ref) != len(alt):
        assert len(alt) == 1 and ref.startswith(alt), "Only the pinned simple deletions are supported"
        start, ref, alt = start + 1, ref[1:], ""
    else:
        assert len(ref) == len(alt) == 1
    offsets = sorted(transcript_offset(reference["exons"], reference["strand"], p)
                     for p in range(start, start + len(ref)))
    assert offsets == list(range(offsets[0], offsets[0] + len(ref)))
    offset = offsets[0]
    if reference["strand"] == "-":
        ref, alt = reverse_complement(ref), reverse_complement(alt)
    cdna = reference["cdna"]
    assert cdna[offset:offset + len(ref)] == ref, (record["gene"], offset, ref)
    mutant_cdna = cdna[:offset] + alt + cdna[offset + len(ref):]
    mutant_protein, has_stop = translate(mutant_cdna[reference["cds_start"]:])
    return dict(reference, variant_offset=offset, oriented_ref=ref, oriented_alt=alt,
                mutant_cdna=mutant_cdna, mutant_protein=mutant_protein,
                mutant_has_stop=has_stop, frameshift=(len(ref) - len(alt)) % 3 != 0)


def aligned_window_start(translation, expected):
    """Independently verify the first RNA base against pinned exon coordinates.

    Do not equate RNA length with reference length across an observed indel.
    No production phase-transfer or genomic-to-transcript helper is used.
    """
    orf = translation.variant_orf
    start = expected["variant_offset"] - len(orf.reference_cdna_sequence_before_variant)
    assert orf.reference_cdna_sequence_before_variant == expected["cdna"][start:expected["variant_offset"]]
    anchors = set()
    reverse = expected["strand"] == "-"
    for read in translation.untrimmed_variant_sequence.reads:
        # RNA coordinate of the ORF's first retained base in this original read.
        query = (len(read.suffix) - orf.variant_cdna_interval_start if reverse
                 else len(read.prefix) - orf.variant_cdna_interval_start)
        if reverse:
            query = len(read.sequence) - 1 - query
        for q0, q1, g0, _ in read.reference_blocks:
            if q0 <= query < q1:
                anchors.add(transcript_offset(expected["exons"], expected["strand"], g0 + query - q0 + 1))
    if any(read.reference_blocks for read in translation.untrimmed_variant_sequence.reads):
        assert anchors == {start}, "RNA/reference left-anchor mismatch"
    return start


def check_translation(translation, expected, protein_length):
    """Check frame, RNA translation, mutation interval and reference-only peptide.

    The last comparison is intentionally separate: real reads can contain
    additional substitutions/indels, so a correctly translated RNA haplotype
    need not equal the reference transcript with only the nominated edit.
    """
    orf = translation.variant_orf
    prefix = orf.variant_cdna_interval_start
    window_start = aligned_window_start(translation, expected)
    assert window_start >= expected["cds_start"]  # all six selected windows are coding
    frame = (expected["cds_start"] - window_start) % 3
    assert orf.offset_to_first_complete_codon == frame
    assert orf.cdna_sequence[prefix:orf.variant_cdna_interval_end] == expected["oriented_alt"]
    assert translation.frameshift == expected["frameshift"]
    actual_aa, actual_stop = translate(orf.cdna_sequence[frame:])
    expected_cdna = expected["mutant_cdna"][window_start:window_start + len(orf.cdna_sequence)]
    expected_aa, expected_stop = translate(expected_cdna[frame:])
    if len(actual_aa) > protein_length:
        actual_aa, actual_stop = actual_aa[:protein_length], False
    if len(expected_aa) > protein_length:
        expected_aa, expected_stop = expected_aa[:protein_length], False
    assert (translation.amino_acids, translation.ends_with_stop_codon) == (actual_aa, actual_stop)
    # Independently locate codons touching the edit; a frame-preserving deletion
    # on a codon boundary marks the zero-width novel peptide junction.
    first_codon = (window_start + frame - expected["cds_start"]) // 3
    local_start = (prefix - frame) // 3
    if expected["frameshift"]:
        local_end = len(actual_aa)
    else:
        local_end = (prefix + len(expected["oriented_alt"]) - frame + 2) // 3
    assert translation.mutation_start_idx == local_start
    assert translation.mutation_end_idx == min(local_end, protein_length)
    assert translation.contains_mutation == (0 < local_start < len(actual_aa)
                                            if local_start == local_end else len(actual_aa) > local_start)
    return {"expected_amino_acids": expected_aa, "expected_stop": expected_stop,
            "matches_expected": (actual_aa, actual_stop) == (expected_aa, expected_stop),
            "protein_start_1based": first_codon + 1}


def reference_genome(cache_directory, directory=REFERENCE):
    """Real PyEnsembl parser/index, restricted to these six pinned transcripts."""
    from pyensembl import Genome

    directory = Path(directory)
    manifest = json.loads((directory / "protein_reference_manifest.json").read_text())
    checksums = {name: metadata["subset_sha256"] for name, metadata in manifest["files"].items()}
    dataset_identity = reference_dataset_identity(directory, checksums)
    genome = Genome(
        reference_name=dataset_identity,
        annotation_name="osteosarc-ensembl-subset", annotation_version=manifest["ensembl_release"],
        gtf_path_or_url=str(directory / "reference.gtf.gz"),
        transcript_fasta_paths_or_urls=[str(directory / "reference.cdna.fa.gz")],
        protein_fasta_paths_or_urls=[str(directory / "reference.pep.fa.gz")],
        copy_local_files_to_cache=True, cache_directory_path=str(Path(cache_directory) / dataset_identity))
    genome.index()
    return genome
