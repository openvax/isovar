# Osteosarc RNA benchmark tools

These commands run from the repository root against the verified offline
osteosarc cache described in
[`tests/data/osteosarc/expansion/README.md`](../../tests/data/osteosarc/expansion/README.md).
They are development tools, outside the installed `isovar` package.

```sh
python -m benchmarks.osteosarc.benchmark --help
python -m benchmarks.osteosarc.collection_benchmark --help
python -m benchmarks.osteosarc.performance_report --help
python -m benchmarks.osteosarc.collection_report --help
python -m benchmarks.osteosarc.coordinate_scan --help
```

The [archived reports and raw measurements](../../docs/archive/experiments/osteosarc-rna-performance/README.md)
retain the original script paths and hashes. Moving the tools does not alter
those measurements or make their original hashes match these files.
