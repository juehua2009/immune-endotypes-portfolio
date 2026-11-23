#!/bin/bash
# Kinchen et al. 2018 (Cell) — GEO GSE114374 supplementary expression matrices.
# Verified on the NCBI FTP listing, Nov 22, 2026.
set -e
DEST=~/Documents/immune-endotypes/data/kinchen_2018
mkdir -p "$DEST"
cd "$DEST"
BASE=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE114nnn/GSE114374/suppl
curl -O $BASE/GSE114374_Human_HC_expression_matrix.txt.gz
curl -O $BASE/GSE114374_Human_UC_expression_matrix.txt.gz
curl -O $BASE/GSE114374_Mouse_DSS_expression_matrix.txt.gz
curl -O $BASE/GSE114374_Mouse_HC_expression_matrix.txt.gz
ls -lh
