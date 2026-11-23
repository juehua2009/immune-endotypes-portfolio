#!/bin/bash
# Martin et al. 2019 (Cell) — GEO GSE134809 supplementary RAW tar (per-sample files).
# Verified on the NCBI FTP listing, Nov 22, 2026: GSE134809_RAW.tar, 579M.
set -e
DEST=~/Documents/immune-endotypes/data/martin_2019
mkdir -p "$DEST"
cd "$DEST"
curl -O https://ftp.ncbi.nlm.nih.gov/geo/series/GSE134nnn/GSE134809/suppl/GSE134809_RAW.tar
tar -tf GSE134809_RAW.tar | head -50
