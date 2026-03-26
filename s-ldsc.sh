#!/bin/bash

LDSC_PATH="path/to/ldsc"
OUT_PATH=""
BASELINE="$LDSC_PATH/baseline_data"
CELL_PATH="superclusters"
GWAS_PATH=""

# Construct annotation path string
ANNOT_PATH="$BASELINE/baselineLD.,$(awk '!/cellnames\.txt/ {print "'$CELL_PATH'/"$0"/baseline."}' cellnames.txt | paste -sd, -)"

FILE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" filenames.txt)

GWAS="${GWAS_PATH}/${FILE}"

TRAIT=$(echo "$FILE" | sed -E 's/\.sumstats\.gz$//')

python2 "$LDSC_PATH/ldsc.py" \
    --h2 "$GWAS" \
    --ref-ld-chr "$ANNOT_PATH" \
    --w-ld-chr "$LDSC_PATH/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC." \
    --overlap-annot \
    --frqfile-chr "$LDSC_PATH/1000G_Phase3_frq/1000G.EUR.QC." \
    --print-coefficients \
    --out "$OUT_PATH/$TRAIT"

echo "Partitioned Heritability analysis for phenotype $FILE completed."
