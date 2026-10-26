"""Read every retained SV3 path from the canonical, manifest-driven audit.

The former one-path calculations now live in isovar.sv_rna_assessment.
The old v1 JSON is a historical snapshot; this consumer does not recompute or
select a highest-support branch.
"""

from .build_unified_sv_manifest import OUTPUT, SV3
from .build_unified_sv_report import snapshots, render


def build(results=OUTPUT / 'results'):
    return [snapshot for _, snapshot in snapshots(results, SV3)]


def main():
    output = OUTPUT / 'SV3.md'
    output.write_text(render(event_id=SV3))
    print(output)


if __name__ == '__main__':
    main()
