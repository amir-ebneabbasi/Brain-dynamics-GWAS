#!/bin/bash

soft_dir=""
out_dir=""
gwas_dir=""

# file format: GWAS_i.txt.fastGWA
trait="state${SLURM_ARRAY_TASK_ID}"
file="${gwas_dir}/GWAS_${SLURM_ARRAY_TASK_ID}.txt.fastGWA"

echo "Running MAGMA for $trait using file $file"

# prepare files
pval_file="${out_dir}/${trait}_pvalue.txt"
N=$(awk 'NR==2 {print $6}' "$file")
awk 'NR > 1 {print $2 "\t" $10}' "$file" > "$pval_file"

for annot_path in "$soft_dir"/hmagma/*.genes.annot; do
    annot_name=$(basename "$annot_path" .genes.annot)
    
    echo "Running MAGMA for $trait with annotation $annot_name"
    
    "$soft_dir"/magma/magma --bfile "$soft_dir/magma/g1000_eur" \
        --pval "$pval_file" N=$N \
        --gene-annot "$soft_dir/hmagma/${annot_name}.genes.annot" \
        --out "$out_dir/h_gene_${annot_name}_${trait}"

    echo "Finished MAGMA for $trait with annotation $annot_name"
done

rm $pval_file
