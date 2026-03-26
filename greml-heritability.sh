#!/bin/bash

GCTA_BIN="gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1"

GRM_MATRIX="unrelated_grm"
QCOVAR="Cov_continous.txt"             
COVAR="Cov_discrete.txt"                           
OUTPUT_DIR=""                 
PHENO_FILE="HMM_for_GWAS.txt"               

# Set the phenotype index from SLURM_ARRAY_TASK_ID
if [ "$SLURM_ARRAY_TASK_ID" -eq 0 ]; then
    # Run without --mpheno
    OUTFILE="${OUTPUT_DIR}/heritability_0.txt"
    $GCTA_BIN --reml --reml-est-fix --reml-pred-rand --grm-cutoff 0.05 \
        --grm "$GRM_MATRIX" \
        --pheno "$PHENO_FILE" \
        --covar "$COVAR" \
        --qcovar "$QCOVAR" \
        --out "$OUTFILE"
else
    # Use --mpheno with the incremented index for other tasks
    MPHENO=$(($SLURM_ARRAY_TASK_ID + 1))
    OUTFILE="${OUTPUT_DIR}/heritability_${SLURM_ARRAY_TASK_ID}.txt"
    $GCTA_BIN --reml --reml-est-fix --reml-pred-rand --grm-cutoff 0.05 \
        --grm "$GRM_MATRIX" \
        --pheno "$PHENO_FILE" \
        --mpheno "$MPHENO" \
        --covar "$COVAR" \
        --qcovar "$QCOVAR" \
        --out "$OUTFILE"
fi

echo "Heritability analysis for phenotype ${MPHENO} completed. Output saved to ${OUTFILE}."
