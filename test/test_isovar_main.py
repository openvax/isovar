from isovar import run_isovar, isovar_results_to_dataframe

def test_isovar_main_to_dataframe():
    results = run_isovar(
        variants="data/b16.f10/b16.vcf",
        alignment_file="data/b16.f10/b16.combined.sorted.bam")
    df = isovar_results_to_dataframe(results)
    print(df)