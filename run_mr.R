#!/usr/bin/env Rscript

run_mr <- function(exp_name, exp_file, out_name, out_file, clump_dir, out_dir) {
	
  cat("Running MR for", exp_name, "->", out_name, "\n")
  
  # =====================================================
  # 1. Load required packages
  # =====================================================
  pkgs <- c("TwoSampleMR", "dplyr", "devtools", "MRPRESSO", "mr.raps")
  
#  for (pkg in pkgs) {
#    if (!require(pkg, character.only = TRUE)) {
#      if (pkg == "TwoSampleMR") {
#        install.packages(pkg, repos = c("https://mrcieu.r-universe.dev","https://cloud.r-project.org"))
#      } else if (pkg == "MRPRESSO") {
#        devtools::install_github("rondolab/MR-PRESSO")
#      } else if (pkg == "mr.raps") {
#        devtools::install_github("qingyuanzhao/mr.raps")
#      } else {
#        install.packages(pkg, repos = "https://cloud.r-project.org")
#      }
#      library(pkg, character.only = TRUE)
#    }
#  }
  
  lapply(pkgs, library, character.only = TRUE)
  
  # =============================
  # 2. Read exposure data
  # =============================
  exp_dat <- read_exposure_data(
    filename = exp_file,
    sep = " ",
    snp_col = "SNP",
    beta_col = "b",
    se_col = "se",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    eaf_col = "freq",
    pval_col = "p",
    samplesize_col = "n"
  )
  exp_dat$exposure <- exp_name
  
  # =============================
  # 3. Filter clumped SNPs
  # =============================
  exp_clumped <- read.table(file.path(clump_dir, paste0("Clump_", exp_name, ".clumped")), header = TRUE)

  num_cols <- ncol(exp_clumped)
  cat("Number of columns in clumped exposure data: ", num_cols, "\n")

  exp_dat <- exp_dat %>% filter(SNP %in% exp_clumped$SNP)
  
  cat("Here is exp_dat after clumping:\n")
  print(head(exp_dat))

  # =============================
  # 4. Read outcome data
  # =============================
  out_dat <- read_outcome_data(
    snps = exp_dat$SNP,
    filename = out_file,
    sep = " ",
    snp_col = "SNP",
    beta_col = "b",
    se_col = "se",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    eaf_col = "freq",
    pval_col = "p",
    samplesize_col = "n"
  )
  out_dat$outcome <- out_name
  
  cat("Here is out_dat:\n")
  print(head(out_dat))

  # =============================
  # 5. Harmonise and  Steiger filtering
  # =============================
  dat <- harmonise_data(exp_dat, out_dat)
  
  # Initial directionality check
  initial_dir_check <- directionality_test(dat)
  cat("Here is initial_dir_check:\n")
  print(initial_dir_check)

  # Ensure you get a proper TRUE/FALSE
  causal_true <- as.logical(initial_dir_check$correct_causal_direction[1])

  # Apply Steiger filtering only if causal direction is FALSE
  if (!causal_true) {
    dat_with_wrong_snps <- steiger_filtering(dat)
    dat <- dat_with_wrong_snps[dat_with_wrong_snps$steiger_dir == TRUE, ]
    removed_steiger <- nrow(dat_with_wrong_snps) - nrow(dat)
  } else {
    removed_steiger <- NA
  }
  
  if (nrow(dat) < 5) stop("Fewer than 5 SNPs remaining after clumping!")

  cat("Number of SNPs removed by Steiger filtering:\n")
  print(removed_steiger)

  cat("Here is harmonised data:\n")
  print(head(dat))

  # =============================
  # 6. Two-sample MR
  # =============================
  mr_results <- mr(dat, method_list=c("mr_ivw", "mr_weighted_median", "mr_egger_regression"))

  cat("Here is mr_results:\n")
  print(mr_results)
  
  # =============================
  # 7. MR-PRESSO
  # =============================
  
  # Try running MR-PRESSO, and retry with 10000 if it fails
  tryCatch({
    presso_results <- mr_presso(
      BetaOutcome     = "beta.outcome",
      BetaExposure    = "beta.exposure",
      SdOutcome       = "se.outcome",
      SdExposure      = "se.exposure",
      OUTLIERtest     = TRUE,
      DISTORTIONtest  = TRUE,
      data            = dat,
      NbDistribution  = 1000,
      SignifThreshold = 0.05
    )
  }, error = function(e) {
    cat("Retrying with 10000...\n")
    presso_results <<- mr_presso(
      BetaOutcome     = "beta.outcome",
      BetaExposure    = "beta.exposure",
      SdOutcome       = "se.outcome",
      SdExposure      = "se.exposure",
      OUTLIERtest     = TRUE,
      DISTORTIONtest  = TRUE,
      data            = dat,
      NbDistribution  = 10000,
      SignifThreshold = 0.05
    )
  })

  # Print the results
  cat("Here is presso_results:\n")
  print(presso_results)

  # =============================
  # 8. Combine MR
  # =============================
  presso_main <- with(presso_results[["Main MR results"]],
                      data.frame(
                        id.exposure = unique(mr_results$id.exposure),
                        id.outcome  = unique(mr_results$id.outcome),
                        exposure    = unique(mr_results$exposure),
                        outcome     = unique(mr_results$outcome),
                        method      = ifelse(`MR Analysis` == "Raw", "MR-PRESSO",
                                             ifelse(`MR Analysis` == "Outlier-corrected", "MR-PRESSO outlier-corrected",
                                                    `MR Analysis`)),
                        nsnp        = unique(mr_results$nsnp),
                        b           = `Causal Estimate`,
                        se          = `Sd`,
                        pval        = `P-value`
                      )
  )
  
  combined_results <- rbind(
    mr_results[, c("id.exposure", "id.outcome", "exposure", "outcome", "method", "nsnp", "b", "se", "pval")],
    presso_main
  )
  
  combined_results$nsnp_removed_steiger <- removed_steiger
  # =============================
  # 9. Sensitivity analyses
  # =============================
  het_results <- mr_heterogeneity(dat)
  cat("Here are het_results:\n")
  print(het_results)
  
  pleio_results <- mr_pleiotropy_test(dat)
  cat("Here are pleio_results:\n")
  print(pleio_results)
  
  direction_results <- directionality_test(dat)
  cat("Here are direction_results:\n")
  print(direction_results)
  
  # =============================
  # 10. Combine sensitivity
  # =============================
  combined_sensitivity <- data.frame(
    exposure = rep(het_results$exposure[1], 5),
    outcome = rep(het_results$outcome[1], 5),
    statistic = c("Q-Egger", "Q-IVW", "Egger intercept", "Global-PRESSO", "Steiger"),
    SE = c(NA, NA, pleio_results$se, NA, NA),
    value = c(
      het_results$Q[1],
      het_results$Q[2],
      pleio_results$egger_intercept,
      presso_results[["MR-PRESSO results"]][["Global Test"]]$RSSobs,
      direction_results$correct_causal_direction
    ),
    p = c(
      het_results$Q_pval[1],
      het_results$Q_pval[2],
      pleio_results$pval,
      presso_results[["MR-PRESSO results"]][["Global Test"]]$Pvalue,
      direction_results$steiger_pval
    )
  )
  
  # =============================
  # 11. Save results
  # =============================
  out_prefix <- file.path(out_dir, paste0(exp_name, "_", out_name))
  
  outputs <- list(
    MR           = combined_results,
    Sensitivity  = combined_sensitivity,
    Harmonized   = dat
  )
  
  for (name in names(outputs)) {
    if (!is.null(outputs[[name]])) {
      write.table(
        outputs[[name]],
        paste0(out_prefix, "_", name, ".tsv"),
        sep = "\t", row.names = FALSE, quote = FALSE
      )
    }
  }
  
  message("✅ Finished successfully for ", exp_name, " → ", out_name)
}
