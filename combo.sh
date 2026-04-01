#!/bin/bash

GCTA="gcta/gcta-1.94.1"

out_dir=""
gwas_dir=""
bfile="fileforclumping"
gene_list="${out_dir}/glist_ensgid_hg19_v40.txt"

trait="state${SLURM_ARRAY_TASK_ID}"
file="${gwas_dir}/GWAS_${SLURM_ARRAY_TASK_ID}.txt.fastGWA"
formatted_file="${out_dir}/GWAS_${SLURM_ARRAY_TASK_ID}_formatted.txt"

echo "Reformatting ${file}..."
awk 'NR==1 {print "SNP\tA1\tA2\tfreq\tBETA\tSE\tP\tN"; next} \
     {print $2"\t"$4"\t"$5"\t"$7"\t"$8"\t"$9"\t"$10"\t"$6}' "${file}" > "${formatted_file}"

echo "Reformatted file saved to ${formatted_file}"

out_prefix="${out_dir}/${trait}"

"${GCTA}" \
  --bfile "${bfile}" \
  --mBAT-combo "${formatted_file}" \
  --mBAT-gene-list "${gene_list}" \
  --mBAT-print-all-p \
  --out "${out_prefix}" \
  --thread-num 10

echo "Finished mBAT-combo for ${formatted_file}"
