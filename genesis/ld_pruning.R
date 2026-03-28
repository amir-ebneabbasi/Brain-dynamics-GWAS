
# -----------------------------
# Load Libraries
# -----------------------------
library(SNPRelate)

# -----------------------------
# Paths
# -----------------------------
gdsfile <- "path/to/genotype_autosomes.gds"
output_dir <- "path/to/output"
pruned_file <- file.path(output_dir, "pruned_snps.txt")

# -----------------------------
# Load GDS file
# -----------------------------
cat("Opening GDS file for LD pruning...\n")
gds <- snpgdsOpen(gdsfile)

# -----------------------------
# Perform LD pruning
# -----------------------------
snpset <- snpgdsLDpruning(gds,
                          method = "corr",
                          slide.max.bp = 10e6,
                          ld.threshold = sqrt(0.1),
                          verbose = T,
                          num.thread = 8)
pruned <- unlist(snpset, use.names = FALSE)

cat("Number of LD-pruned SNPs:", length(pruned), "\n")

# -----------------------------
# Save pruned SNP IDs
# -----------------------------
write.table(pruned, file = pruned_file, row.names = FALSE, col.names = FALSE, quote = TRUE)
cat("LD-pruned SNP list saved to:", pruned_file, "\n")

snpgdsClose(gds)
cat("LD pruning complete.\n")
