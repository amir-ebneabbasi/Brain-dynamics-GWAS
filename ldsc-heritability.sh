#!/bin/bash

LDSC_PATH=""
LOC_PATH=""

python2 "$LDSC_PATH/munge_sumstats.py" \
    --sumstats "$LOC_PATH/GWAS/GWAS_${SLURM_ARRAY_TASK_ID}.txt.fastGWA" \
    --out "$LOC_PATH/munge/M_${SLURM_ARRAY_TASK_ID}" \
    --merge-alleles "$LDSC_PATH/w_hm3.snplist" \
    --chunksize 500000 

echo "Munge analysis for phenotype ${SLURM_ARRAY_TASK_ID} completed"


python2 "$LDSC_PATH/ldsc.py" \
    --h2 "$LOC_PATH/munge/M_${SLURM_ARRAY_TASK_ID}.sumstats.gz" \
    --ref-ld-chr "$LDSC_PATH/eur_w_ld_chr/" \
    --w-ld-chr "$LDSC_PATH/eur_w_ld_chr/" \
    --out "$LOC_PATH/ldcs_heritability/H_${SLURM_ARRAY_TASK_ID}" \

echo "Heritability analysis for phenotype ${SLURM_ARRAY_TASK_ID} completed"
