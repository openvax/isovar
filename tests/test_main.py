from isovar import run_isovar, isovar_results_to_dataframe
from nose.tools import eq_
from testing_helpers import data_path

def test_isovar_main_to_dataframe():
    results = run_isovar(
        variants=data_path("data/b16.f10/b16.vcf"),
        alignment_file=data_path("data/b16.f10/b16.combined.sorted.bam"))
    df = isovar_results_to_dataframe(results)
    print(df)
    eq_(len(df), 4)
    # B16 test data has 2/4 variants with enough coverage
    # to translate protein sequences
    eq_(df["passes_all_filters"].sum(), 2)

