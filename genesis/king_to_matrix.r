# -----------------------------
# Load Libraries
# -----------------------------
library(GWASTools)
library(SNPRelate)
library(GENESIS)

# -----------------------------
# Paths
# -----------------------------
dir <- "path/to/dir"
gdsfile <- file.path(dir, "genotype_autosomes.gds")
king_file <- file.path(dir, "eur_king.kin0")
geno <- GdsGenotypeReader(filename = gdsfile)
genoData <- GenotypeData(geno)

# -----------------------------
# Get Ids
# -----------------------------
iids <- as.character(getScanID(genoData))

# -----------------------------
# Convert King output to matrix
# -----------------------------
KINGmat <- kingToMatrix(
  king_file,
  sample.include = iids,
  thresh = 2^(-11/2),
  estimator = "Kinship"
)

# -----------------------------
# Save the matrix
# -----------------------------
saveRDS(KINGmat, file = file.path(dir, "eur_kinship_matrix.rds"))
cat("Kinship matrix saved as eur_kinship_matrix.rds\n")
