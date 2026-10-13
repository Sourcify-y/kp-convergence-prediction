#!/usr/bin/env bash
set -euo pipefail
find accessions_8318.txt data scripts figures requirements.txt software_versions.txt -type f | sort | xargs sha256sum > SHA256SUMS
