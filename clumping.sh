#!/bin/bash

module load plink/1.9
echo "PLINK 1.9 loaded."

GWAS_PATH=""
FILE_FOR_CLUMPING="fileforclumping"
CLUMP_OUTPUT=""

plink \
  --bfile "$FILE_FOR_CLUMPING" \
  --clump "$GWAS_PATH/GWAS_${SLURM_ARRAY_TASK_ID}.txt.fastGWA" \
  --clump-field P \
  --clump-p1 5e-8 \
  --clump-p2 1 \
  --clump-r2 0.1 \
  --clump-kb 1000 \
  --threads 15 \
  --out "$CLUMP_OUTPUT/Clump_${SLURM_ARRAY_TASK_ID}"

echo "PLINK clumping completed for task ID: ${SLURM_ARRAY_TASK_ID}"
