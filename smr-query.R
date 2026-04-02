#!/usr/bin/env python3

import csv
import subprocess
from pathlib import Path

base = "path/to/base"
clumped_csv = "path/to/all_clumped.csv"
SMR_Tool = "path/to/smr"

eqtl_loc = f"{base}/smr/eqtl/BrainMeta_cis_eqtl_summary"
sqtl_loc = f"{base}/smr/sqtl/BrainMeta_cis_sqtl_summary"
mqtl_loc = f"{base}/smr/mqtl/Brain-mMeta/Brain-mMeta"
caqtl_loc = f"{base}/smr/caqtl/Bryois_caQTL_summary/bryois_NatCommun_2018_50kb_cQTLs"

out_dir = Path(f"{base}/smr/query_results")
out_dir.mkdir(parents=True, exist_ok=True)


def run_smr(qtl_file: Path, qtl_name: str, snp: str, gwas_id: str):
    besd_path = Path(str(qtl_file) + ".besd")

    if besd_path.exists():
        print(f"Running SMR on {qtl_name}: {qtl_file}")
        out_path = out_dir / f"GWAS_{gwas_id}_{snp}_{qtl_name}"

        cmd = [
            SMR_Tool,
            "--beqtl-summary", str(qtl_file),
            "--query", "5.0e-8",
            "--snp", snp,
            "--out", str(out_path)
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode == 0:
            print(f"{qtl_name} run successfully for SNP {snp}")
        else:
            print(f"{qtl_name} failed for SNP {snp}")
            if result.stdout.strip():
                print(result.stdout.strip())
            if result.stderr.strip():
                print(result.stderr.strip())
    else:
        print(f"Missing {qtl_name} file: {qtl_file} (skipped)")


# run function
with open(clumped_csv, newline='') as f:
    reader = csv.DictReader(f)
    for row in reader:
        chr_ = row["CHR"]
        snp = row["SNP"]
        gwas_id = row["GWAS_id"].strip()

        print("==========================================================")
        print(f"Processing SNP: {snp} (chr{chr_}) for GWAS: {gwas_id}")
        print("==========================================================")

        eqtl_file = Path(f"{eqtl_loc}/BrainMeta_cis_eQTL_chr{chr_}")
        sqtl_file = Path(f"{sqtl_loc}/BrainMeta_cis_sQTL_chr{chr_}")
        mqtl_file = Path(mqtl_loc)
        caqtl_file = Path(caqtl_loc)

        run_smr(eqtl_file, "eQTL", snp, gwas_id)
        run_smr(sqtl_file, "sQTL", snp, gwas_id)
        run_smr(mqtl_file, "mQTL", snp, gwas_id)
        run_smr(caqtl_file, "caQTL", snp, gwas_id)
