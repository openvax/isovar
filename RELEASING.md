# Releasing Isovar

1. On a feature branch, bump `__version__` in `isovar/__init__.py` as part of the PR
   ([semver](https://semver.org/): minor for new features or behavior changes,
   patch for fixes and docs). Add an entry to [CHANGELOG.md](CHANGELOG.md) for any
   change that can alter results or break callers.
2. Merge the PR once CI passes. CI requires the version to be newer than the base
   branch's and not yet released; if another PR released that version first, bump again.
3. From a clean, up-to-date `master`, run `./deploy.sh`. It runs `./lint.sh` and
   `./test.sh`, builds the sdist and wheel, uploads them to PyPI, and tags the release.
   It refuses to run from another branch or with uncommitted changes.
