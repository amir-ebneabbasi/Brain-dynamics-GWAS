
# -----------------------------
# Load Libraries
# -----------------------------
library(GWASTools)
library(SNPRelate)
library(BiocParallel)
library(GENESIS)

# -----------------------------
# Paths
# -----------------------------
# OLD
output_dir <- "path/to/output"
gdsfile <- file.path(output_dir, "genotype_autosomes.gds")
pruned_file <- file.path(output_dir, "pruned_snps.txt")
king_file <- file.path(output_dir, "kinship_matrix.rds")

# NEW
mypc_file <- file.path(output_dir, "mypc.rds")

# -----------------------------
# Open GDS
# -----------------------------
gds <- snpgdsOpen(gdsfile)

# -----------------------------
# Get All Sample IDs
# -----------------------------
all_ids <- read.gdsn(index.gdsn(gds, "sample.id"))

# -----------------------------
# Load pruned SNPs
# -----------------------------
pruned <- scan(pruned_file, what = character())

# -----------------------------
# Load KING Matrix
# -----------------------------
KINGmat <- readRDS(king_file)

# -----------------------------
# Sample Partitioning
# -----------------------------
#sampset <- pcairPartition(
#    kinobj = KINGmat,
#    kin.thresh = 2^(-11/2),
#    divobj = KINGmat,
#    div.thresh = -2^(-11/2)
#)

# -----------------------------
# Save sampset
# -----------------------------
#saveRDS(sampset, file = sampset_file, compress = TRUE)
#cat("sampset saved to:", sampset_file, "\n")

# -----------------------------
# Open GDS for PC analysis
# -----------------------------
snpgdsClose(gds)
gds_reader <- GdsGenotypeReader(filename = gdsfile)
genoData <- GenotypeData(gds_reader)

# -----------------------------
# Run PC-Air
# -----------------------------
mypc <- pcair(genoData, kinobj = KINGmat, divobj = KINGmat, 
          snp.include = pruned, num.cores=8)

# -----------------------------
# Save PC Results
# -----------------------------
saveRDS(mypc, file = mypc_file, compress = TRUE)
cat("PC Results completed and saved to:", mypc_file, "\n")
