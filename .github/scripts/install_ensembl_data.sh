#!/usr/bin/env bash
set -euo pipefail

install_release() {
  local release="$1"
  local species="$2"
  local mirror="$3"
  local attempt
  local delay

  for attempt in 1 2 3; do
    if pyensembl install \
        --release "${release}" \
        --species "${species}" \
        --custom-mirror "${mirror}"; then
      return
    fi
    if [[ "${attempt}" -eq 3 ]]; then
      echo "Failed to install Ensembl release ${release} after ${attempt} attempts" >&2
      return 1
    fi
    delay=$((attempt * 15))
    echo "Ensembl release ${release} install attempt ${attempt} failed; retrying in ${delay}s" >&2
    sleep "${delay}"
  done
}

install_release 87 homo_sapiens \
  https://github.com/openvax/ensembl-data/releases/download/GRCh38.87/
install_release 75 human \
  https://github.com/openvax/ensembl-data/releases/download/GRCh37.75/
install_release 102 mouse \
  https://github.com/openvax/ensembl-data/releases/download/GRCm38.102/
