#!/bin/bash

conda activate path/to/polyfun
echo "Polyfun environment activated: $CONDA_PREFIX"

BASE=""
FINEMAP="$BASE/finemapping"
PARQUET="$FINEMAP/parquet"
SNPVAR="$FINEMAP/snpvar_files"
OUTPUT="$FINEMAP/output"
GWAS="$BASE/GWAS"

mkdir -p "$PARQUET" "$SNPVAR" "$OUTPUT"

ID="GWAS_${SLURM_ARRAY_TASK_ID}"

# Munge summary statistics
python "$FINEMAP/polyfun/munge_polyfun_sumstats.py" \
  --sumstats "$GWAS/${ID}.txt.fastGWA" \
  --out "$PARQUET/${ID}.parquet" \
  --min-info 0.6 \
  --min-maf 0.001

echo "Munge analysis for phenotype ${SLURM_ARRAY_TASK_ID} completed"

# Extract SNP variance
python "$FINEMAP/polyfun/extract_snpvar.py" \
  --sumstats "$PARQUET/${ID}.parquet" \
  --allow-missing \
  --out "$SNPVAR/${ID}_snpvar.parquet"

echo "snpvar analysis for phenotype ${SLURM_ARRAY_TASK_ID} completed"

# Create finemapper jobs
python "$FINEMAP/polyfun/create_finemapper_jobs.py" \
  --geno "path/to/geno" \
  --n 52000 \
  --sumstats "$SNPVAR/${ID}_snpvar.parquet" \
  --method susie \
  --pvalue-cutoff 5e-8 \
  --max-num-causal 5 \
  --out-prefix "$OUTPUT/${ID}" \
  --jobs-file "$OUTPUT/${ID}_jobs.txt"

echo "finemapper jobs for phenotype ${SLURM_ARRAY_TASK_ID} completed"
