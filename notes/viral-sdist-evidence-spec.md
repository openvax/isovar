# Package the reproducible antigen/viral validation evidence

The runtime viral command passed its installed-wheel experiment, but the source
archive omits its offline EBV fixture and the recent antigen/viral benchmark
recipes. Add only the relevant small fixture, scripts, validation notes and
checksum records to MANIFEST.in. Bump to 1.56.1; runtime APIs and algorithms are
unchanged. Build both artifacts, inspect the sdist, unpack it and run the offline
viral tests using its own fixture. Then run lint/full tests/CI, merge and deploy
with the standard script. Historical scientific evidence remains pinned to the
version/run it describes.
