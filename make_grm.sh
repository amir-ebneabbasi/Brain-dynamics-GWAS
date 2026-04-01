#!/bin/bash

BASE="path/to/base"
CHR_PATH="$BASE/chromosome_path.txt"
SNPLIST="$BASE/snplist.txt"
OUTDIR="$BASE/output"

GCTA="path/to/gcta-1.95.0-linux-kernel-3-x86_64/gcta64"

$GCTA \
--mbfile $CHR_PATH \
--extract $SNPLIST \
--make-grm \
--sparse-cutoff 0.05 \
--thread-num 2 \
--out $OUTDIR/sp0.05_autosome_grm

echo "Job completed successfully: GRM saved to $OUTDIR/sp0.05_autosome_grm"

