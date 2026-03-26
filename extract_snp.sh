#!/bin/bash

module load qctool
module load plink/1.9

IMPUTED_DIR=""
MFI_DIR=""
OUT_DIR=""
SNPFILE="$OUT_DIR/snpfile.txt"

mkdir -p "$OUT_DIR"

# Loop through each SNP
while read -r SNP; do
    echo "-----------------------------------"
    echo "Processing SNP: $SNP"

    MFI_FILE=$(grep -l -w "$SNP" "$MFI_DIR"/ukb_mfi_chr*_v3.txt | head -n1)

    if [[ -z "$MFI_FILE" ]]; then
        echo "SNP $SNP not found in any chromosome. Skipping."
        continue
    fi

    CHR=$(basename "$MFI_FILE" | sed -E 's/ukb_mfi_chr([0-9]+)_v3.txt/\1/')
    echo "SNP $SNP found on chromosome $CHR"

    BGEN_FILE="$IMPUTED_DIR/ukb_imp_chr${CHR}_v3.bgen"
    SAMPLE_FILE="$IMPUTED_DIR/ukb20904_imp_chr${CHR}_v3_s487334.sample"
    OUT_PREFIX="$OUT_DIR/extract-${SNP}"

    qctool -g "$BGEN_FILE" -s "$SAMPLE_FILE" -og "${OUT_PREFIX}.bed" -incl-rsids <(echo "$SNP")
    if [[ $? -ne 0 ]]; then
        echo "Failed to extract $SNP from chromosome $CHR"
        continue
    fi

    echo "Extraction completed for $SNP"

    plink --bfile "$OUT_PREFIX" --recodeA --out "$OUT_PREFIX"
    if [[ $? -ne 0 ]]; then
        echo "PLINK recoding failed for $SNP"
        continue
    fi

    echo "Finished SNP: $SNP"
done < "$SNPFILE"

