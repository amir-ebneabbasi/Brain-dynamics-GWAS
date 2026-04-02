#!/usr/bin/env Rscript

library(readr)

run_smr <- function(qtl_file, qtl_name, snp, gwas_id,
                    smr_tool, out_dir) {
  
  besd_path <- paste0(qtl_file, ".besd")
  
  if (file.exists(besd_path)) {
    cat(sprintf("Running SMR on %s: %s\n", qtl_name, qtl_file))
    
    out_path <- file.path(out_dir, paste0("GWAS_", gwas_id, "_", snp, "_", qtl_name))
    
    cmd <- c(
      "--beqtl-summary", qtl_file,
      "--query", "5.0e-8",
      "--snp", snp,
      "--out", out_path
    )
    
    result <- system2(smr_tool, args = cmd, stdout = TRUE, stderr = TRUE)
    exit_status <- attr(result, "status")
    
    if (is.null(exit_status) || exit_status == 0) {
      cat(sprintf("%s run successfully for SNP %s\n", qtl_name, snp))
    } else {
      cat(sprintf("%s failed for SNP %s\n", qtl_name, snp))
      if (length(result) > 0) {
        cat(paste(result, collapse = "\n"), "\n")
      }
    }
    
  } else {
    cat(sprintf("Missing %s file: %s (skipped)\n", qtl_name, qtl_file))
  }
}


# # how to run

# base <- "path/to/base"
# smr_tool <- "path/to/smr"
# clumped_csv <- "path/to/all_clumped.csv"
# eqtl_loc <- file.path(base, "smr/eqtl/BrainMeta_cis_eqtl_summary")
# sqtl_loc <- file.path(base, "smr/sqtl/BrainMeta_cis_sqtl_summary")
# mqtl_loc <- file.path(base, "smr/mqtl/Brain-mMeta/Brain-mMeta")
# caqtl_loc <- file.path(base, "smr/caqtl/Bryois_caQTL_summary/bryois_NatCommun_2018_50kb_cQTLs")
# out_dir <- file.path(base, "smr/query_results")
# dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


# df <- read_csv(clumped_csv, show_col_types = FALSE)

# for (i in seq_len(nrow(df))) {
#   row <- df[i, ]
  
#   chr_ <- row$CHR
#   snp <- row$SNP
#   gwas_id <- trimws(row$GWAS_id)
  
#   cat("==========================================================\n")
#   cat(sprintf("Processing SNP: %s (chr%s) for GWAS: %s\n", snp, chr_, gwas_id))
#   cat("==========================================================\n")
  
#   eqtl_file <- file.path(eqtl_loc, paste0("BrainMeta_cis_eQTL_chr", chr_))
#   sqtl_file <- file.path(sqtl_loc, paste0("BrainMeta_cis_sQTL_chr", chr_))
#   mqtl_file <- mqtl_loc
#   caqtl_file <- caqtl_loc
  
#   run_smr(eqtl_file, "eQTL", snp, gwas_id, smr_tool, out_dir)
#   run_smr(sqtl_file, "sQTL", snp, gwas_id, smr_tool, out_dir)
#   run_smr(mqtl_file, "mQTL", snp, gwas_id, smr_tool, out_dir)
#   run_smr(caqtl_file, "caQTL", snp, gwas_id, smr_tool, out_dir)
# }
