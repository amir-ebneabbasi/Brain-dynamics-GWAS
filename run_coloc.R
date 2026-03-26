#!/usr/bin/env Rscript

run_coloc <- function(line_num, gcta_address, ld_reference, exp_name, exp_file, out_name, out_file, out_dir, ld_snps) {
  
  cat("\n=======================================================================\n")
  cat(" Line", line_num, ": Running coloc for:", exp_name, "→", out_name, "\n")
  cat("=======================================================================\n\n")
  
  library(ieugwasr)
  library(TwoSampleMR)
  library(coloc)
  library(data.table)
  library(dplyr)
  library(devtools)
  
  # ===========================================
  # 2. Read exposure data
  # ===========================================
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

  # ===========================================
  # 3. Filter based on provided LD SNPs
  # ===========================================
  cat("→ Filtering exp_dat by provided LD SNPs...\n")
  
  ld_snp_list <- unlist(strsplit(ld_snps, ","))
  ld_snp_list <- trimws(ld_snp_list)
  
  exp_dat_filtered <- exp_dat %>% filter(SNP %in% ld_snp_list)
  
  cat("Number of LD SNPs:", length(ld_snp_list), "\n")
  cat("Number of rows in exp_dat_filtered:", nrow(exp_dat_filtered), "\n\n")
  
  if (nrow(exp_dat_filtered) == 0) {
   stop("No matching SNPs found in exp_dat_filtered")
  }

  # ===========================================
  # 4. Read outcome data
  # ===========================================
  out_dat <- read_outcome_data(
    snps = exp_dat_filtered$SNP,
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
  cat("→ Outcome data loaded with", nrow(out_dat), "variants.\n\n")
  
  # ===========================================
  # 5. Harmonise
  # ===========================================
  dat <- harmonise_data(exp_dat_filtered, out_dat) %>% na.omit()
  if (nrow(dat) == 0) stop("No SNPs left after harmonisation.")
  cat("→ Harmonised data size:", nrow(dat), "\n")

  # ===========================================
  # 6. Calculate LD matrix (to be aligned with effect allele of data)
  # ===========================================
  snp_info <- dat %>%
    dplyr::select(SNP, effect_allele.exposure)
  
  snp_file <- file.path(out_dir, paste0("temp_snps_", line_num, "_", exp_name, "_", out_name, ".txt"))
  
  write.table(snp_info,
              file = snp_file,
              quote = FALSE,
              sep = "\t",
              row.names = FALSE,
              col.names = FALSE)
  
  cat("Saved SNP and effect allele list to:", snp_file, "\n")
  
  cat("Running GCTA to update reference alleles ...\n")
  
  ld_prefix <- file.path(out_dir, paste0("temp_ld_", line_num, "_", exp_name, "_", out_name))
  
  gcta_cmd <- paste(
    shQuote(gcta_address),
    "--bfile", shQuote(ld_reference),
    "--extract", shQuote(snp_file),
    "--update-ref-allele", shQuote(snp_file),
    "--recode",
    "--out", shQuote(ld_prefix)
  )
  
  cat("Executing:", gcta_cmd, "\n")
  system(gcta_cmd)
  
  ld_file <- paste0(ld_prefix, ".xmat.gz")
  ld_update_log <- paste0(ld_prefix, ".log")
  
  if (!file.exists(ld_file)) {
    stop("LD file not found after running GCTA: ", ld_file)
  }

  snp_header <- scan(ld_file, what = "", nlines = 1, quiet = TRUE)
  snp_header <- snp_header[!snp_header %in% c("FID", "IID")]
  
  cat("Total SNPs in updated LD file:", length(snp_header), "\n")
  cat("First 5 SNPs:\n")
  print(head(snp_header, 5))
  
  filtered_dat <- dat[dat$SNP %in% snp_header, ]
  cat("\nNumber of dat rows after filtering by header of updated LD:", nrow(filtered_dat), "\n")
  cat("→ First five SNPs:\n")
  print(head(filtered_dat$SNP, 5))
  
  cat("\nReading genotype matrix ...\n")
  
  geno_matrix <- read.table(ld_file,
                            skip = 2,
                            stringsAsFactors = FALSE)
  
  # Remove first two columns (FID IID)
  geno_matrix <- geno_matrix[, -c(1,2)]
  
  colnames(geno_matrix) <- snp_header
  
  cat("\nComputing LD matrix (pairwise correlation)...\n")
  ld_matrix <- cor(as.matrix(geno_matrix), use = "pairwise.complete.obs")
  
  # Replace NA with 0 (no variation → undefined LD; prevents coloc.susie errors)
  ld_matrix[is.na(ld_matrix)] <- 0

  # Reorder SNPs to match coloc input (data and LD must have identical order)
  snp_order <- filtered_dat$SNP
  ld_matrix <- ld_matrix[snp_order, snp_order]
  
  cat("\nLD matrix dimensions:", dim(ld_matrix), "\n")
  cat("Preview of LD matrix (top-left 5x5):\n")
  print(ld_matrix[1:min(5, nrow(ld_matrix)), 1:min(5, ncol(ld_matrix))])

  # ===========================================
  # 7. Coloc analysis
  # ===========================================
  cat("→ Running coloc analysis...\n")
  
  # Get sample size
  get_mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  
  N_exposure <- get_mode(filtered_dat$samplesize.exposure)
  N_outcome  <- get_mode(filtered_dat$samplesize.outcome)
  
  cat("  Sample size for exposure: ", N_exposure, "\n")
  cat("  Sample size for outcome: ", N_outcome, "\n")
  
  # Get variable type
  infer_type <- function(name) {
    cc_keywords <- c("ADHD","Alcohol","ALS","ALZ","Anxiety","ASD","Bipolar",
                     "DLB","Insomnia","MDD","OCD","Panic","PD","PTSD","Schizophrenia",
                     "Stroke","Substance","FTD","VasDem")
    if (any(grepl(paste(cc_keywords, collapse = "|"), name, ignore.case = TRUE))) {
      return("cc")
    } else {
      return("quant")
    }
  }
  
  type_exposure <- infer_type(exp_name)
  type_outcome  <- infer_type(out_name)
  
  cat("  exposure type: ", type_exposure, "\n")
  cat("  outcome type: ", type_outcome, "\n")
  
  # Prepare datasets
  dataset_exposure <- list(
    snp     = filtered_dat$SNP,
    beta    = filtered_dat$beta.exposure,
    varbeta = filtered_dat$se.exposure^2,
    N       = N_exposure,
    MAF     = filtered_dat$eaf.exposure,
    LD      = ld_matrix,
    type    = type_exposure
  )
  
  dataset_outcome <- list(
    snp     = filtered_dat$SNP,
    beta    = filtered_dat$beta.outcome,
    varbeta = filtered_dat$se.outcome^2,
    N       = N_outcome,
    MAF     = filtered_dat$eaf.outcome,
    LD      = ld_matrix,
    type    = type_outcome
  )
  
  # LD check
  cat("→ Checking LD alignment in exposure...\n")
  check_ld <- tryCatch({
    check_alignment(dataset_exposure, thr = 0.2, do_plot = FALSE)
  }, error = function(e) {
    cat("⚠ LD alignment check failed with error:\n", e$message, "\n")
    NULL
  })
  
  if (!is.null(check_ld)) {
    cat("Correlation of LD and beta values in exposure is:\n")
    print(check_ld)
  } else {
    cat("⚠ Correlation of LD and beta values cannot be reported.\n")
  }
  
  # Prefix for saving
  out_prefix <- file.path(out_dir, paste0("results_", line_num, "_", exp_name, "_", out_name))
  
  # ===========================================
  # Run coloc.abf
  # ===========================================
  cat("→ Checking datasets for coloc.abf...\n")
  check_exp <- check_dataset(dataset_exposure)
  check_out <- check_dataset(dataset_outcome)
  
  if (is.null(check_exp) && is.null(check_out)) {
    cat("✓ Datasets OK for coloc.abf. Running analysis...\n")
    res_abf <- coloc.abf(dataset1 = dataset_exposure, dataset2 = dataset_outcome)
    
    if (!is.null(res_abf)) {
      # Compute credible set
      o <- order(res_abf$results$SNP.PP.H4, decreasing = TRUE)
      cs <- cumsum(res_abf$results$SNP.PP.H4[o])
      w <- which(cs > 0.95)[1]
      credible_set <- res_abf$results[o, ][1:w, ]$snp
      
      res_abf_df <- as.data.frame(as.list(res_abf[[1]]))
      res_abf_df$credible <- paste(credible_set, collapse = ",")
      res_abf_df$Var1 <- exp_name
      res_abf_df$Var2 <- out_name
      
      write.table(res_abf_df,
                  paste0(out_prefix, "_abf.txt"),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      cat("✓ coloc.abf results saved.\n")
    }
  } else {
    cat("✗ Dataset check failed for coloc.abf. Issues:\n")
    if (!is.null(check_abf_exp)) print(check_abf_exp)
    if (!is.null(check_abf_out)) print(check_abf_out)
  }
  
  # ===========================================
  # Run coloc.susie
  # ===========================================
  if (!is.null(check_ld)) {
    
    if (is.null(check_exp) && is.null(check_out)) {
      cat("✓ Datasets OK for coloc.susie. Running analysis...\n")
      res_susie <- coloc.susie(dataset1 = dataset_exposure, dataset2 = dataset_outcome)
      
      if (!is.null(res_susie$summary)) {
        res_susie_df <- as.data.frame(res_susie$summary)
        res_susie_df$Var1 <- exp_name
        res_susie_df$Var2 <- out_name
        res_susie_df$LD_Beta_corr_in_Var1 <- check_ld
        
        write.table(res_susie_df,
                    paste0(out_prefix, "_susie.txt"),
                    sep = "\t", quote = FALSE, row.names = FALSE)
        cat("✓ coloc.susie results saved.\n")
      }
    } else {
      cat("✗ Dataset check failed for coloc.susie. Issues:\n")
      if (!is.null(check_susie_exp)) print(check_susie_exp)
      if (!is.null(check_susie_out)) print(check_susie_out)
    }
  } else {
    cat("⚠ Skipping coloc.susie due to LD alignment check failure.\n")
  }
  
  # ===========================================
  # Clean up temporary files
  # ===========================================
  if (file.exists(snp_file)) file.remove(snp_file)
  if (file.exists(ld_file)) file.remove(ld_file)
  if (file.exists(ld_update_log)) file.remove(ld_update_log)
  
  cat("\n✅ Finished successfully for", exp_name, "→", out_name, "\n")
}
