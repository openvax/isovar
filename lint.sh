#!/bin/bash
set -o errexit

ruff check isovar/ tests/

echo 'Passes ruff check'
