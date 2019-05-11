from isovar.isovar import Isovar

def test_isovar_create_dataframe():
    isovar_instance = Isovar()
    df = isovar_instance.create_dataframe(
        variants="data/b16.f10/b16.vcf",
        alignment_file="data/b16.f10/b16.combined.sorted.bam")
    print(df)