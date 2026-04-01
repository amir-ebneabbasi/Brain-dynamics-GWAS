run_hyprcoloc_pipeline <- function(input_dir, out_dir, p_threshold = 5e-8) {

  # ===========================================
  # 1. Load libraries
  # ===========================================
  library(data.table)
  library(dplyr)
  library(mapgen)
  library(ieugwasr)
  library(tidyr)
  library(stringr)
  library(gtools)
  library(hyprcoloc)
  library(readr)
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # ===========================================
  # 2. Loading LD blocks
  # ===========================================
  message("Loading LD blocks...")
  LD_Blocks <- readRDS(system.file('extdata', 'LD.blocks.EUR.hg19.rds', package = 'mapgen')) %>%
    mutate(
      chr = as.integer(chr),
      start = as.integer(start),
      end = as.integer(end)
    )

  # ===========================================
  # 3. Extracting significant SNPs
  # ===========================================
  message("Extracting significant SNPs...")
  
  all_sig_snps <- data.frame()
  files <- list.files(input_dir, pattern = "\\.fastGWA$", full.names = TRUE)
  
  trait_names <- tools::file_path_sans_ext(basename(files))
  trait_names <- gsub("\\.fastGWA$", "", trait_names)
  
  for (filepath in files) {
    filename <- basename(filepath)
    message("Processing ", filename)
    
    df <- fread(filepath, data.table = FALSE) %>%
      filter(!is.na(CHR), !is.na(POS), !is.na(P)) %>%
      mutate(
        CHR = as.integer(CHR),
        POS = as.integer(POS)
      )
    
    sig_snps <- df %>% filter(P <= p_threshold)
    all_sig_snps <- bind_rows(all_sig_snps, sig_snps)
  }
  
  all_sig_snps <- all_sig_snps %>% distinct(SNP, .keep_all = TRUE)

  # ===========================================
  # 4. Mapping SNPs to LD blocks
  # ===========================================
  message("Mapping SNPs to LD blocks...")
  
  nearby_blocks <- data.frame()
  
  for (i in seq_len(nrow(all_sig_snps))) {
    row <- all_sig_snps[i, ]
    block <- LD_Blocks %>%
      filter(chr == row$CHR,
             start <= row$POS,
             end >= row$POS)
    
    nearby_blocks <- bind_rows(nearby_blocks, block)
  }
  
  nearby_blocks <- nearby_blocks %>% distinct(locus, .keep_all = TRUE)
  
  write_tsv(nearby_blocks, file.path(out_dir, "nearby_block.txt"))

  # ===========================================
  # 5. Building beta/SE matrices
  # ===========================================
  message("Building beta/SE matrices...")
  
  for (block_idx in seq_len(nrow(nearby_blocks))) {
    
    block <- nearby_blocks[block_idx, ]
    message("Processing block ", block_idx, " (locus ", block$locus, ")")
    
    beta_block <- NULL
    se_block <- NULL
    
    for (filepath in files) {
      
      df <- read_tsv(filepath, show_col_types = FALSE)
      
      if (ncol(df) == 1) {
        df <- read_delim(filepath, delim = " ", show_col_types = FALSE)
      }
      
      file_name <- tools::file_path_sans_ext(basename(filepath))
      trait_name <- gsub("\\.fastGWA$", "", file_name)
      
      block_df <- df %>%
        filter(CHR == block$chr,
               POS >= block$start,
               POS <= block$end)
      
      if (nrow(block_df) == 0) next
      
      beta_col <- paste0(trait_name, "_b")
      se_col   <- paste0(trait_name, "_se")
      
      beta_tmp <- block_df %>%
        select(SNP, BETA) %>%
        rename(!!beta_col := BETA) %>%
        mutate(!!beta_col := ifelse(!!sym(beta_col) == 0, 1e-8, !!sym(beta_col)))
      
      se_tmp <- block_df %>%
        select(SNP, SE) %>%
        rename(!!se_col := SE) %>%
        mutate(!!se_col := ifelse(!!sym(se_col) == 0, 1e-8, !!sym(se_col)))
      
      beta_block <- if (is.null(beta_block)) beta_tmp else full_join(beta_block, beta_tmp, by = "SNP")
      se_block   <- if (is.null(se_block)) se_tmp else full_join(se_block, se_tmp, by = "SNP")
    }
    
    if (!is.null(beta_block)) {
      write_tsv(beta_block,
                file.path(out_dir, paste0("beta_block_", block_idx, "_locus_", block$locus, ".txt")))
    }
    
    if (!is.null(se_block)) {
      write_tsv(se_block,
                file.path(out_dir, paste0("se_block_", block_idx, "_locus_", block$locus, ".txt")))
    }
  }

  # ===========================================
  # 6. Running hyprcoloc
  # ===========================================
  message("Running hyprcoloc...")
  
  beta_files <- list.files(out_dir, pattern = "^beta_block_.*\\.txt$", full.names = TRUE)
  
  for (beta_file in beta_files) {
    
    message("Processing ", basename(beta_file))
    
    se_file <- gsub("beta_", "se_", beta_file)
    
    if (!file.exists(se_file)) {
      message("Missing SE file for ", beta_file)
      next
    }
    
    beta_df <- read.table(beta_file, header = TRUE)
    se_df   <- read.table(se_file, header = TRUE)
    
    if (!all(beta_df$SNP == se_df$SNP)) {
      stop("SNP mismatch in ", beta_file)
    }
    
    rsid <- beta_df$SNP
    
    rownames(beta_df) <- rsid
    rownames(se_df)   <- rsid
    
    beta_mat <- as.matrix(beta_df[, setdiff(colnames(beta_df), "SNP")])
    se_mat   <- as.matrix(se_df[, setdiff(colnames(se_df), "SNP")])
    
    traits <- gsub("_b$", "", colnames(beta_mat))
    
    tryCatch({
      res <- hyprcoloc(beta_mat, se_mat,
                       trait.names = traits,
                       snp.id = rsid,
                       uniform.priors = FALSE)
      
      write.csv(res[["results"]],
                file = file.path(out_dir,
                                 gsub("^beta_", "results_",
                                      gsub("\\.txt$", ".csv", basename(beta_file)))),
                row.names = FALSE)
      
    }, error = function(e) {
      message("Error in ", beta_file, ": ", e$message)
    })
  }
  
  message("Pipeline completed successfully.")
}
