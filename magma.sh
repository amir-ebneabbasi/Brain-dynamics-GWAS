#!/bin/bash

soft_dir=""
out_dir=""
gwas_dir=""

# file format: GWAS_i.txt.fastGWA
trait="state${SLURM_ARRAY_TASK_ID}"
file="${gwas_dir}/GWAS_${SLURM_ARRAY_TASK_ID}.txt.fastGWA"

echo "Running MAGMA for $trait using file $file"

# prepare files
loc_file="${out_dir}/${trait}_loc.txt"
pval_file="${out_dir}/${trait}_pvalue.txt"
awk 'NR > 1 {print $2 "\t" $1 "\t" $3}' "$file" > "$loc_file"

# Annotate genes
"$soft_dir"/magma --annotate window=35,10 \
    --snp-loc "$loc_file" \
    --gene-loc "$out_dir/NCBI37.3.gene.loc.extendedMHCexcluded" \
    --out "$out_dir/annot_${trait}"
rm "$loc_file"

# Gene analysis
awk 'NR > 1 {print $2 "\t" $10}' "$file" > "$pval_file"
N=$(awk 'NR==2 {print $6}' "$file")
"$soft_dir"/magma --bfile "$soft_dir/g1000_eur" \
    --pval "$pval_file" N=$N \
    --gene-annot "$out_dir/annot_${trait}.genes.annot" \
    --out "$out_dir/gene_${trait}"
rm "$pval_file"

# Gene set enrichment
"$soft_dir"/magma --gene-results "$out_dir/gene_${trait}.genes.raw" \
    --model direction=greater \
    --gene-covar "$out_dir/Siletti_l2_conti_specificity_matrix.txt" \
    --out "$out_dir/enrich_${trait}"

echo "Finished MAGMA for $trait"
