#!/usr/bin/env bash
set -euo pipefail

output_dir="${1:-work/public_cache/predictors}"
mkdir -p "$output_dir"

curl -fL --retry 3 --continue-at - \
  -o "$output_dir/revel_G.bw" \
  "https://hgdownload.soe.ucsc.edu/gbdb/hg38/revel/g.bw"

curl -fL --retry 3 --continue-at - \
  -o "$output_dir/alphaMissense_G.bw" \
  "https://hgdownload.soe.ucsc.edu/gbdb/hg38/alphaMissense/g.bw"
