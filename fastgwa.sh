#!/bin/bash

GCTA="gcta/gcta-1.94.1"

GRM_MATRIX="sp0.05_autosome_grm"
MBFILE="GWAS_mbfile.txt"
QCOVAR="Cov_continous.txt"             
COVAR="Cov_discrete.txt"                           
OUTPUT_DIR=""                 
PHENO_FILE="HMM_for_GWAS.txt"               

# Set the phenotype index from SLURM_ARRAY_TASK_ID
if [ "$SLURM_ARRAY_TASK_ID" -eq 0 ]; then
    # Run without --mpheno
    OUTFILE="${OUTPUT_DIR}/GWAS_0.txt"
    $GCTA --fastGWA-mlm \
        --grm-sparse "$GRM_MATRIX" \
	    --mbfile "$MBFILE" \
        --pheno "$PHENO_FILE" \
        --covar "$COVAR" \
        --qcovar "$QCOVAR" \
        --out "$OUTFILE"
else
    # Use --mpheno with the incremented index for other tasks
    MPHENO=$(($SLURM_ARRAY_TASK_ID + 1))
    OUTFILE="${OUTPUT_DIR}/GWAS_${SLURM_ARRAY_TASK_ID}.txt"
    $GCTA --fastGWA-mlm \
        --grm-sparse "$GRM_MATRIX" \
	    --mbfile "$MBFILE" \
        --pheno "$PHENO_FILE" \
        --mpheno "$MPHENO" \
        --covar "$COVAR" \
        --qcovar "$QCOVAR" \
        --out "$OUTFILE"
fi

echo "GWAS analysis for phenotype ${MPHENO} completed. Output saved to ${OUTFILE}."
