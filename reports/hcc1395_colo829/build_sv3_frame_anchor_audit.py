"""Thin consumer of the shared all-path SV3 assessment.

The historical v1 frame-anchor JSON documents the previous selected-path audit.
New evidence and frame decisions are owned by the library and canonical
manifest snapshots, including every retained alternative.
"""

from .build_sv3_path_fusion import build, main

__all__ = ['build', 'main']

if __name__ == '__main__':
    main()
