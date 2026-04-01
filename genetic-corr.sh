#!/bin/bash

LDSC_PATH=""
STATE_PATH=""
DISEASE_PATH=""
OUT_PATH=""

mkdir -p "$OUT_PATH"

PAIR=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" rg_pairs.txt)
STATE_FILE=$(echo "$PAIR" | cut -d' ' -f1)
DISEASE_FILE=$(echo "$PAIR" | cut -d' ' -f2)

STATE_GWAS="$STATE_PATH/${STATE_FILE}.sumstats.gz"
DISEASE_GWAS="$DISEASE_PATH/${DISEASE_FILE}.sumstats.gz"

echo "Running RG for: $STATE_FILE vs $DISEASE_FILE"

python2 "$LDSC_PATH/ldsc.py" \
    --rg "$STATE_GWAS","$DISEASE_GWAS" \
    --ref-ld-chr "$LDSC_PATH/eur_w_ld_chr/" \
    --w-ld-chr "$LDSC_PATH/eur_w_ld_chr/" \
    --out "$OUT_PATH/${STATE_FILE}_${DISEASE_FILE}"