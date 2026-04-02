run_gene_to_gwas_pipeline <- function(input_csv, output_csv, delay = 0.3) {
  
  genesymbol2gwas <- function(gene_input) {
    url <- paste0(
      "https://www.ebi.ac.uk/gwas/api/search/downloads?",
      "q=ensemblMappedGenes:", gene_input,
      "&pvalfilter=&orfilter=&betafilter=&datefilter=&genomicfilter=",
      "&genotypingfilter[]=&traitfilter[]=&dateaddedfilter=",
      "&facet=association&efo=true"
    )
    
    df <- tryCatch({
      read.delim(url, sep = "\t", stringsAsFactors = FALSE)
    }, error = function(e) {
      message(sprintf("Error reading %s: %s", gene_input, e$message))
      return(NULL)
    })
    
    return(df)
  }
  
  # Load all genes
  df <- read.csv(input_csv, stringsAsFactors = FALSE)
  
  results <- list()
  counter <- 1
  
  for (i in seq_len(nrow(df))) {
    
    gene_input <- df$Gene[i]
    gwas_id_input <- df$GWAS_ID[i]
    
    d <- genesymbol2gwas(gene_input)
    
    if (!is.null(d) && nrow(d) > 0) {
      
      names(d) <- gsub("[./ ]", "_", names(d))
      
      cols_needed <- c("GWAS_ID", "GENE", "DISEASE_TRAIT", "STUDY_ACCESSION")
      cols_present <- intersect(cols_needed, names(d))
      
      if (length(cols_present) == 0) {
        message(sprintf("No expected columns found for %s", gene_input))
        next
      }
      
      d <- d[, cols_present, drop = FALSE]
      
      # Add metadata
      d$gene_input <- gene_input
      d$gwas_id_input <- gwas_id_input
      
      results[[counter]] <- d
      counter <- counter + 1
      
      message(sprintf("Done: %s (%d rows)", gene_input, nrow(d)))
    } else {
      message(sprintf("No data for %s", gene_input))
    }
    
    Sys.sleep(delay)  # avoid API throttling
  }
  
  # Combine results
  df_gwas <- if (length(results) > 0) do.call(rbind, results) else data.frame()
  
  df_gwas <- unique(df_gwas)
  
  write.csv(df_gwas, output_csv, row.names = FALSE)
  
  message(sprintf("Final shape: %d rows, %d cols", nrow(df_gwas), ncol(df_gwas)))
  message(sprintf("Saved to: %s", output_csv))
  
  return(df_gwas)
}