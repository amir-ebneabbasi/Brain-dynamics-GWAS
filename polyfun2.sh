#!/bin/bash

conda activate path/to/polyfun
echo "Polyfun environment activated: $CONDA_PREFIX"

job_list_file="All_GWAS_jobs.txt"
geno_base=""

line=$((SLURM_ARRAY_TASK_ID + 1))
job_command=$(sed -n "${line}p" "$job_list_file")

if [[ -n "$job_command" ]]; then
  # Extract chromosome number from --chr argument
  chr=$(echo "$job_command" | grep -oP '(?<=--chr )\d+')

  # Replace the --geno argument with the correct path
  job_command=$(echo "$job_command" | sed -E "s|--geno [^ ]+|--geno ${geno_base}${chr}|")

  echo "Running modified command: $job_command"
  eval "$job_command"
else
  echo "No command found at line $line"
fi