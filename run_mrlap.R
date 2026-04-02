#!/usr/bin/env Rscript

# =====================================================
# TwoSampleMR + LDSC + MRlap correction
# =====================================================

run_mr <- function(exp_name,
                   exp_file,
                   out_name,
                   out_file,
                   clump_dir,
                   out_dir,
                   MR_threshold = NULL) {

  cat("Running MR for", exp_name, "→", out_name, "\n")

  # =============================
  # 1. Load packages
  # =============================
  pkgs <- c("TwoSampleMR", "dplyr", "data.table")
  suppressMessages(lapply(pkgs, library, character.only = TRUE))

  # =============================
  # 2. Read full GWAS summary stats (for LDSC / MRlap)
  # =============================
  cat("Reading full exposure data", "\n")
  exp_full <- fread(exp_file)

  cat("Reading full outcome data", "\n")
  out_full <- fread(out_file)

  required_cols <- c("SNP","CHR","POS","BETA","SE","A1","A2","N","P")

  missing_exp <- setdiff(required_cols, colnames(exp_full))
  missing_out <- setdiff(required_cols, colnames(out_full))

  if (length(missing_exp) > 0) {
    cat("Missing in exp_full:", paste(missing_exp, collapse = ", "), "\n")
  }

  if (length(missing_out) > 0) {
    cat("Missing in out_full:", paste(missing_out, collapse = ", "), "\n")
  }

  stopifnot(length(missing_exp) == 0)
  stopifnot(length(missing_out) == 0)

  # =============================
  # 3. Read exposure and outcome for MR
  # =============================
  exp_dat <- read_exposure_data(
    filename = exp_file,
    snp_col = "SNP",
    chr_col = "CHR",
    pos_col = "POS",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    eaf_col = "MAF",
    pval_col = "P",
    samplesize_col = "N"
  )
  exp_dat$exposure <- exp_name

  out_dat <- read_outcome_data(
    filename = out_file,
    snp_col = "SNP",
    chr_col = "CHR",
    pos_col = "POS",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    eaf_col = "MAF",
    pval_col = "P",
    samplesize_col = "N"
  )
  out_dat$outcome <- out_name

  # =============================
  # 4. Filter on clumped SNPs
  # =============================
  clump_file <- file.path(clump_dir, paste0("Clump_", exp_name, ".clumped"))
  exp_clumped <- fread(clump_file)
  iv_snps <- exp_clumped$SNP

  exp_dat <- exp_dat %>% filter(SNP %in% iv_snps)
  out_dat <- out_dat %>% filter(SNP %in% iv_snps)
  cat("Number of IVs after clumping:", nrow(exp_dat), "\n")

  # =============================
  # 5. Harmonise + Steiger filtering
  # =============================
  dat <- harmonise_data(exp_dat, out_dat)
  dir_test <- directionality_test(dat)
  causal_ok <- as.logical(dir_test$correct_causal_direction[1])

  if (!causal_ok) {
    dat0 <- dat
    dat <- steiger_filtering(dat)
    dat <- dat[dat$steiger_dir == TRUE, ]
    cat("Steiger filtering removed:", nrow(dat0) - nrow(dat), "SNPs\n")
  }

  if (nrow(dat) < 5)
    stop("Fewer than 5 SNPs remain after filtering")

  # =============================
  # 6. Two-sample MR (IVW)
  # =============================
  mr_results <- mr(dat, method_list = "mr_ivw")
  ivw_row <- mr_results[mr_results$method == "Inverse variance weighted", ]

  # =============================
  # 7. Source MRlap functions
  # =============================
  source('https://github.com/n-mounier/MRlap/raw/refs/heads/master/R/get_correction.R')
  source('https://github.com/n-mounier/MRlap/raw/refs/heads/master/R/run_LDSC.R')

  # =============================
  # 8. Run LDSC on full GWAS
  # =============================
  ldsc_params <- run_LDSC(
    exposure_data = exp_full,
    exposure_name = exp_name,
    outcome_data  = out_full,
    outcome_name  = out_name,
    ld            = "path/to/ldsc/eur_w_ld_chr",
    hm3           = "path/to/ldsc/eur_w_ld_chr/w_hm3.snplist",
    save_logfiles = FALSE,
    verbose       = TRUE
  )

  # =============================
  # 9. Determine MR threshold
  # =============================
  if (is.null(MR_threshold)) {
    MR_threshold <- if (grepl("state", exp_name, ignore.case = TRUE)) 5e-6 else 5e-8
  }

  # =============================
  # 10. Run MRlap correction
  # =============================
  mrlap_res <- get_correction(
    IVs = dat %>%
      select(beta.exposure, se.exposure) %>%
      mutate(across(everything(), as.numeric)) %>%
      rename(std_beta.exp = beta.exposure,
           std_SE.exp   = se.exposure),

    lambda         = as.numeric(ldsc_params$lambda)[1],
    lambda_se      = as.numeric(ldsc_params$lambda_se)[1],
    h2_LDSC        = as.numeric(ldsc_params$h2_LDSC)[1],
    h2_LDSC_se     = as.numeric(ldsc_params$h2_LDSC_se)[1],

    alpha_obs      = ivw_row$b,
    alpha_obs_se   = ivw_row$se,

    n_exp          = max(as.numeric(dat$samplesize.exposure), na.rm = TRUE),
    n_out          = max(as.numeric(dat$samplesize.outcome),  na.rm = TRUE),

    MR_threshold   = MR_threshold,
    verbose        = TRUE
  )

  # =============================
  # 11. Save final results
  # =============================
  if (!is.null(mrlap_res)) {
    final_results <- data.frame(
      id.exposure = dat$id.exposure[1],
      id.outcome  = dat$id.outcome[1],
      exposure    = exp_name,
      outcome     = out_name,
      method      = 'MRlap-corrected IVW',
      nsnp        = mr_results$nsnp,
      b           = mrlap_res$alpha_corrected,
      se          = mrlap_res$alpha_corrected_se,
      pval        = 2 * stats::pnorm(-abs(mrlap_res$alpha_corrected / mrlap_res$alpha_corrected_se))
    )

    out_file_name <- file.path(out_dir, paste0(exp_name, "_", out_name, "_MRlap.tsv"))
    fwrite(final_results, out_file_name, sep = "\t")

    message("Result is successfully written for ", exp_name, " → ", out_name)
  } else {
    warning("No results saved")
  }
}
