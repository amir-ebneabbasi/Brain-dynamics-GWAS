#!/bin/bash

LDSC_PATH=""
LOC_PATH=""

for gwas_file in "$LOC_PATH/data/"*.{tsv,txt}; do

    base_name=$(basename "$gwas_file" | sed 's/\.[^.]*$//')

    if python2 "$LDSC_PATH/munge_sumstats.py" \
        --sumstats "$gwas_file" \
        --out "$LOC_PATH/data/munged/$base_name" \
        --merge-alleles "$LDSC_PATH/w_hm3.snplist" \
        --chunksize 500000; then

        echo "Munge analysis for phenotype ${base_name} completed"
    else
        echo "ERROR: Munge failed for ${base_name}" >&2
    fi
done