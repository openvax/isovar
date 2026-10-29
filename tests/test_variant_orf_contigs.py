"""#622: supplied ORF loci use exact names and Varcode's ambiguity rules."""

import pytest
from orfdata import ORFModel, ORFRecord
from pyensembl import Genome
from pyensembl.species import find_species_by_name
from varcode import Variant

from isovar import AlleleRead, ProteinSequenceCreator, ReadEvidence, reconstruct_variant_orfs


@pytest.fixture
def orf_reference(tmp_path):
    genomes = []

    def build(contigs):
        sequence = "C" * 15 + "ATG" + "GCC" * 10 + "TAA" + "G" * 15
        gtf, fasta = tmp_path / "reference.gtf", tmp_path / "reference.fa"
        rows = []
        for index, contig in enumerate(contigs):
            attrs = (f'gene_id "G{index}"; gene_name "TEST{index}"; transcript_id "T{index}"; '
                     'gene_biotype "lncRNA"; transcript_biotype "lncRNA";')
            rows.extend(f'{contig}\ttest\t{kind}\t100\t165\t.\t+\t.\t{attrs}{extra}\n'
                        for kind, extra in [("gene", ""), ("transcript", ""),
                                            ("exon", f' exon_id "E{index}"; exon_number "1";')])
        gtf.write_text("".join(rows))
        fasta.write_text("".join(f">T{index}\n{sequence}\n" for index in range(len(contigs))))
        genome = Genome("test_assembly", "test_annotation", "1", gtf_path_or_url=str(gtf),
                        transcript_fasta_paths_or_urls=[str(fasta)], cache_directory_path=str(tmp_path))
        genome.species = find_species_by_name("homo_sapiens")
        genome.index()
        genomes.append(genome)
        references = [ORFRecord.from_model(ORFModel(
            f"orf{index}", f"T{index}", f"G{index}", f"TEST{index}", contig, "+", ((115, 150),),
            "M" + "A" * 10, "lncRNA", "test_assembly", 9606, "homo_sapiens", source="synthetic"
        )).bind(genome.transcript_by_id(f"T{index}")) for index, contig in enumerate(contigs)]
        return genome, references, sequence

    yield build
    for genome in genomes:
        genome.close()


@pytest.mark.parametrize("contigs,query,status,selected", [
    (("chrUn_MTfoo", "chrUn_Mfoo"), "chrUn_MTfoo", "exact", "T0"),
    (("chrUn_MTfoo", "chrUn_Mfoo"), "chrUn_Mfoo", "exact", "T1"),
    (("MT",), "chrM", "alias", "T0"),
    (("chrM", "MT"), "chrM", "exact", "T0"),
    (("1", "CHR1"), "chr1", "ambiguous", None),
    (("1", "CHR1"), "CHR1", "exact", "T1"),
    (("HPV16",), "HPV16", "exact", "T0"),
    (("hs38d1",), "hs38d1", "exact", "T0"),
    (("KI270442.1",), "chrUn_KI270442v1", "alias", "T0"),
    (("KI270442.1",), "KI270442.2", "missing", None),
    (("chrUn_Mfoo",), "chrUn_MTfoo", "missing", None),
    (("1",), "unavailable", "missing", None),
])
def test_supplied_reference_contigs_are_not_equated_by_string_replacement(
        orf_reference, contigs, query, status, selected):
    genome, references, sequence = orf_reference(contigs)
    variant = Variant(query, 134, sequence[34], "T", genome, normalize_contig_names=False)
    reads = [AlleleRead(sequence[:34], "T", sequence[35:], name=f"read{i}",
                        reference_blocks=((0, len(sequence), 99, 165),)) for i in range(2)]
    evidence = ReadEvidence.from_variant_and_allele_reads(variant, reads)
    result = reconstruct_variant_orfs(
        variant, evidence, references, sample_id="sample", source="library",
        protein_sequence_creator=ProteinSequenceCreator(protein_sequence_preference="support"))
    assert result["reference_contig_resolution"]["status"] == status
    assert result["reference_contig_resolution"]["original_name"] == query
    assert {candidate["transcript_id"] for candidate in result["protein_candidates"]} == (
        set() if selected is None else {selected})
    assert result["rna_candidates"]  # unresolved reference vocabulary does not erase RNA
    if selected is None:
        assert all(model["reason"] == "unresolved_reference_contig" for model in result["models"])
        assert all(row["state"] == "no_supplied_reference_path" for row in result["rna_candidates"])
    else:
        assert all(row["compatible_transcript_ids"] == [selected] for row in result["rna_candidates"])


def test_empty_reference_vocabulary_is_unavailable_not_zero_rna(orf_reference):
    genome, _, sequence = orf_reference(("1",))
    variant = Variant("1", 134, sequence[34], "T", genome)
    evidence = ReadEvidence.from_variant_and_allele_reads(variant, [
        AlleleRead(sequence[:34], "T", sequence[35:], name=f"read{i}") for i in range(2)])
    result = reconstruct_variant_orfs(variant, evidence, [], sample_id="sample", source="library")
    assert result["reference_contig_resolution"]["status"] == "reference_unavailable"
    assert result["assembly_state"] == "assembled"
