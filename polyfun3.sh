#!/bin/bash

conda activate path/to/polyfun
echo "Polyfun environment activated: $CONDA_PREFIX"

BASE=""
FINEMAP="$BASE/finemapping"
SNPVAR="$FINEMAP/snpvar_files"
OUTPUT="$FINEMAP/output"

ID="GWAS_${SLURM_ARRAY_TASK_ID}"

python "$FINEMAP/polyfun/aggregate_finemapper_results.py" \
    --out-prefix "$OUTPUT/${ID}" \
    --sumstats "$SNPVAR/${ID}_snpvar.parquet" \
    --allow-missing-jobs \
    --pvalue-cutoff 5e-8 \
    --out "$OUTPUT/${ID}_polyfun_agg.txt.gz"

echo "Aggregate jobs for phenotype ${SLURM_ARRAY_TASK_ID} completed"