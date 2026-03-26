#!/bin/bash

base=""
smr_loc="$base/smr"
eqtl_loc="$smr_loc/eqtl/BrainMeta_cis_eqtl_summary"
sqtl_loc="$smr_loc/sqtl/BrainMeta_cis_sqtl_summary"
mqtl_loc="$smr_loc/mqtl/Brain-mMeta/Brain-mMeta"
caqtl_loc="$smr_loc/caqtl/Bryois_caQTL_summary/bryois_NatCommun_2018_50kb_cQTLs"
LD_loc="LD_files"
GWAS_ORIG="$base/GWAS/GWAS_${SLURM_ARRAY_TASK_ID}.txt.fastGWA"
GWAS_SMR="$smr_loc/GWAS_${SLURM_ARRAY_TASK_ID}_smr.txt"

awk 'NR==1 {print "SNP\tA1\tA2\tfreq\tb\tse\tp\tn"; next} \
     {print $2"\t"$4"\t"$5"\t"$7"\t"$8"\t"$9"\t"$10"\t"$6}' $GWAS_ORIG > $GWAS_SMR

SMR_Tool="path/to/smr"

# Loop through chromosomes
for i in {1..22}; do
    echo "Running SMR for chr${i} in $GWAS_SMR ..."

    # eQTL
    $SMR_Tool \
        --bfile $LD_loc/ukbchr_v2_r2correct_${i}_forLD_v2 \
        --gwas-summary $GWAS_SMR \
        --beqtl-summary $eqtl_loc/BrainMeta_cis_eQTL_chr${i} \
        --out $smr_loc/eqtl_chr${i}_GWAS_${SLURM_ARRAY_TASK_ID} \
        --thread-num 10
    echo "Done eQTL chr${i} for $GWAS_SMR"

    # sQTL
    $SMR_Tool \
        --bfile $LD_loc/ukbchr_v2_r2correct_${i}_forLD_v2 \
        --gwas-summary $GWAS_SMR \
        --beqtl-summary $sqtl_loc/BrainMeta_cis_sQTL_chr${i} \
        --out $smr_loc/sqtl_chr${i}_GWAS_${SLURM_ARRAY_TASK_ID} \
        --thread-num 10
    echo "Done sQTL chr${i} for $GWAS_SMR"

    # mQTL
    $SMR_Tool \
        --bfile $LD_loc/ukbchr_v2_r2correct_${i}_forLD_v2 \
        --gwas-summary $GWAS_SMR \
        --beqtl-summary $mqtl_loc \
        --out $smr_loc/mqtl_chr${i}_GWAS_${SLURM_ARRAY_TASK_ID} \
        --thread-num 10
    echo "Done mQTL chr${i} for $GWAS_SMR"

    # caQTL
    $SMR_Tool \
        --bfile $LD_loc/ukbchr_v2_r2correct_${i}_forLD_v2 \
        --gwas-summary $GWAS_SMR \
        --beqtl-summary $caqtl_loc \
        --out $smr_loc/caqtl_chr${i}_GWAS_${SLURM_ARRAY_TASK_ID} \
        --thread-num 10
    echo "Done caQTL chr${i} for $GWAS_SMR"
done

rm $GWAS_SMR
