
# -----------------------------
# Load Libraries
# -----------------------------
library(dplyr)
library(tibble)
library(readr)
library(GENESIS)
library(methods)
library(GWASTools)
library(SNPRelate)
library(BiocParallel)
library(SeqVarTools)

# -----------------------------
# Base directory
# -----------------------------
dir <- "path/to/dir"
outdir <- file.path(dir, "genesis", "output")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# -----------------------------
# File paths
# -----------------------------
mydat_file <- file.path(dir, "HMM_for_GENESIS.txt")
grm_file   <- file.path(dir, "genesis/eur", "mypcrelate_sparse.rds")
gds_file   <- file.path(dir, "genesis/eur", "genotype_autosomes.gds")

# -----------------------------
# Check files exist
# -----------------------------
stopifnot(file.exists(mydat_file))
stopifnot(file.exists(grm_file))
stopifnot(file.exists(gds_file))

# -----------------------------
# Make ScanAnnotationDataFrame
# -----------------------------
mydat <- read.table(mydat_file, header = TRUE, stringsAsFactors = FALSE)
scanAnnot <- ScanAnnotationDataFrame(mydat)

# -----------------------------
# SLURM task ID
# -----------------------------
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
outcome_var <- paste0("State_", task_id)
cat("Running for outcome:", outcome_var, "\n")

# -----------------------------
# Ensure numeric covariates
# -----------------------------
# NOTE:
# age2, age_gender, and age2_gender are removed to preclude multicollinearity in the null model
covariate_names <- c(
  "scansite", "gender", "age",
  "euler", "fd", "fd_max",
  paste0("PC_AiR_", 1:16)
)

# Convert covariates to numeric
for (cov in covariate_names) {
  scanAnnot[[cov]] <- as.numeric(as.character(scanAnnot[[cov]]))
}

# Check for NAs after conversion
na_counts <- sapply(covariate_names, function(x) sum(is.na(scanAnnot[[x]])))
if (any(na_counts > 0)) {
  warning(
    "Some covariates have NAs after conversion:\n",
    paste(names(na_counts)[na_counts > 0], collapse = ", ")
  )
}

# -----------------------------
# Fit null mixed model
# -----------------------------
grm <- readRDS(grm_file)

nullmod <- fitNullModel(
  scanAnnot,
  outcome = outcome_var,
  covars  = covariate_names,
  cov.mat = grm,
  family  = "gaussian"
)

cat("Null model fitted for", outcome_var, "\n")

# -----------------------------
# Load genotype data
# -----------------------------
geno <- GdsGenotypeReader(gds_file)
genoData = GenotypeData(geno)

genoiter <- GenotypeBlockIterator(
  genoData,
  snpBlock = 5000
)

# -----------------------------
# Single SNP association testing
# -----------------------------
assoc <- assocTestSingle(
  genoiter,
  null.model = nullmod,
  BPPARAM = MulticoreParam(workers = 16)
)

assoc_rds_file <- file.path(outdir, paste0("assoc_", task_id, ".rds"))
saveRDS(assoc, file = assoc_rds_file)
cat("Association results saved to:", assoc_rds_file, "\n")

# -----------------------------
# Heritability estimation
# -----------------------------
h2 <- varCompCI(nullmod, prop = TRUE)

h2_file <- file.path(outdir, paste0("heritability_", task_id, ".rds"))
saveRDS(h2, file = h2_file)
cat("Heritability estimates saved to:", h2_file, "\n")

# -----------------------------
# Format GWAS results (fastGWA format)
# -----------------------------
assoc_tsv_file <- file.path(outdir, paste0("GWAS_", task_id, ".txt"))

alleles <- tibble(
  SNP = getSnpID(genoData),
  A1  = getAlleleA(genoData),
  A2  = getAlleleB(genoData)
)

assoc %>%
  as_tibble() %>%
  select(-MAC) %>%
  rename(
    SNP  = variant.id,
    CHR  = chr,
    POS  = pos,
    N    = n.obs,
    AF1  = freq,
    BETA = Est,
    SE   = Est.SE,
    P    = Score.pval
  ) %>%
  inner_join(alleles, by = "SNP") %>%
  mutate(N = N / 2) %>%
  select(CHR, SNP, POS, A1, A2, N, AF1, BETA, SE, P) %>%
  write_tsv(assoc_tsv_file)

cat("Formatted GWAS output saved to:", assoc_tsv_file, "\n")

# -----------------------------
# Close GDS connection
# -----------------------------
close(genoData)
cat("GDS connection closed\n")
