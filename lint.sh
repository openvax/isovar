#!/usr/bin/env bash
set -euo pipefail

resolve_python() {
  local candidate
  local resolved_python
  if [[ -n "${PYTHON:-}" ]]; then
    candidate="${PYTHON}"
  elif [[ -n "${VIRTUAL_ENV:-}" && -x "${VIRTUAL_ENV}/bin/python" ]]; then
    candidate="${VIRTUAL_ENV}/bin/python"
  elif [[ -x ".venv/bin/python" ]]; then
    candidate=".venv/bin/python"
  else
    candidate="python3"
  fi

  if ! resolved_python="$(command -v "${candidate}" 2>/dev/null)" || \
      [[ ! -f "${resolved_python}" || ! -x "${resolved_python}" ]]; then
    echo "Python interpreter not found or not executable: ${candidate}" >&2
    exit 1
  fi
  if [[ "${resolved_python}" != /* ]]; then
    resolved_python="${PWD}/${resolved_python#./}"
  fi
  PYTHON="${resolved_python}"
  export PYTHON
}

resolve_python
"${PYTHON}" -m ruff check isovar/ tests/

echo 'Passes ruff check'
