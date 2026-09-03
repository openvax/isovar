#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: ./deploy.sh [--dry-run] [version]" >&2
}

DRY_RUN=0
VERSION_INPUT=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run)
      DRY_RUN=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      if [[ -z "${VERSION_INPUT}" ]]; then
        VERSION_INPUT="$1"
        shift
      else
        echo "Unexpected argument: $1" >&2
        usage
        exit 1
      fi
      ;;
  esac
done

CURRENT_BRANCH="$(git rev-parse --abbrev-ref HEAD)"
if [[ "${CURRENT_BRANCH}" != "main" && "${CURRENT_BRANCH}" != "master" ]]; then
  echo "Deploys are only allowed from main or master (current: ${CURRENT_BRANCH})." >&2
  exit 1
fi

if [[ -n "$(git status --porcelain)" ]]; then
  echo "Working tree not clean; commit or stash changes before deploy." >&2
  exit 1
fi

PYTHON="${PYTHON:-python3}"
if [[ -x ".venv/bin/python" ]]; then
  PYTHON=".venv/bin/python"
fi

CURRENT_VERSION="$(
  "$PYTHON" - <<'PY'
import re
from pathlib import Path

text = Path("isovar/__init__.py").read_text()
match = re.search(r'__version__ = "([^"]+)"', text)
if not match:
    raise SystemExit("Failed to read isovar/__init__.py; unexpected format.")
print(match.group(1))
PY
)"

VERSION="${CURRENT_VERSION}"
if [[ -n "${VERSION_INPUT}" ]]; then
  VERSION="${VERSION_INPUT#v}"
fi
if [[ ! "${CURRENT_VERSION}" =~ ^(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)$ ]]; then
  echo "Invalid current version in isovar/__init__.py: ${CURRENT_VERSION}" >&2
  exit 1
fi
if [[ ! "${VERSION}" =~ ^(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)$ ]]; then
  echo "Invalid release version: ${VERSION}" >&2
  exit 1
fi
if [[ -n "${VERSION_INPUT}" && "${VERSION}" != "${CURRENT_VERSION}" ]]; then
  "$PYTHON" - "${CURRENT_VERSION}" "${VERSION}" <<'PY'
import sys

current = tuple(int(part) for part in sys.argv[1].split("."))
target = tuple(int(part) for part in sys.argv[2].split("."))
if target <= current:
    raise SystemExit(
        "Release version must increase: %s -> %s" % (sys.argv[1], sys.argv[2])
    )
PY
fi

TAG="v${VERSION}"
LOCAL_TAG_EXISTS=0
if git rev-parse -q --verify "refs/tags/${TAG}" >/dev/null; then
  LOCAL_TAG_EXISTS=1
  TAG_COMMIT="$(git rev-parse "${TAG}^{commit}")"
  HEAD_COMMIT="$(git rev-parse HEAD)"
  if [[ "${VERSION}" != "${CURRENT_VERSION}" || \
        "${TAG_COMMIT}" != "${HEAD_COMMIT}" ]]; then
    echo "Tag ${TAG} already exists locally at a different release commit; aborting." >&2
    exit 1
  fi
  echo "Resuming release from existing local tag ${TAG}."
fi

./lint.sh
./test.sh

if [[ "${DRY_RUN}" -eq 1 ]]; then
  if [[ -n "${VERSION_INPUT}" && "${VERSION}" != "${CURRENT_VERSION}" ]]; then
    echo "Dry run: would bump ${CURRENT_VERSION} to ${VERSION}."
  fi
  echo "Dry run: local gates passed for ${TAG}; no files or remote state changed."
  exit 0
fi

set +e
git ls-remote --exit-code --tags origin "refs/tags/${TAG}" >/dev/null 2>&1
REMOTE_TAG_STATUS=$?
set -e
if [[ "${REMOTE_TAG_STATUS}" -eq 0 ]]; then
  echo "Tag ${TAG} already exists on origin; aborting." >&2
  exit 1
elif [[ "${REMOTE_TAG_STATUS}" -ne 2 ]]; then
  echo "Could not determine whether ${TAG} exists on origin; aborting." >&2
  exit 1
fi

if [[ -n "${VERSION_INPUT}" && "${VERSION}" != "${CURRENT_VERSION}" ]]; then
  VERSION="${VERSION}" "$PYTHON" - <<'PY'
import os
import re
from pathlib import Path

version = os.environ["VERSION"]
path = Path("isovar/__init__.py")
text = path.read_text()
new_text, n = re.subn(
    r'__version__ = "([^"]+)"',
    '__version__ = "%s"' % version,
    text,
    count=1,
)
if n != 1:
    raise SystemExit("Failed to update isovar/__init__.py; unexpected format.")
path.write_text(new_text)
PY
  git add isovar/__init__.py
  git commit -m "Bump version to ${VERSION}"
fi

# Always push before publishing. This also recovers cleanly if an earlier
# attempt committed the bump locally but failed during its push.
git push origin "${CURRENT_BRANCH}"

"$PYTHON" -m pip install --upgrade build twine
rm -rf dist build
"$PYTHON" -m build
"$PYTHON" -m twine check dist/*
"$PYTHON" release_upload.py \
  --project isovar \
  --version "${VERSION}" \
  dist/*

if [[ "${LOCAL_TAG_EXISTS}" -eq 0 ]]; then
  git tag "${TAG}"
fi
git push origin "${TAG}"
echo "Deployed isovar ${VERSION}: https://pypi.org/project/isovar/${VERSION}/"
