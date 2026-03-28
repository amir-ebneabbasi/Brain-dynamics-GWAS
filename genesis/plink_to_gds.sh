# load libraries
library(SNPRelate)

# Input and output paths
input_prefix <- "path/to/autosomes"
output_dir <- "path/to/output"

# Construct file names
bed.fn <- sprintf("%s.bed", input_prefix)
bim.fn <- sprintf("%s.bim", input_prefix)
fam.fn <- sprintf("%s.fam", input_prefix)
out.gdsfn <- sprintf("%s/genotype_autosomes.gds", output_dir)

# Check files exist before conversion
for (f in c(bed.fn, bim.fn, fam.fn)) {
  if (!file.exists(f)) stop(paste("❌ File not found:", f))
}

# Convert merged PLINK files to GDS
cat("Converting merged autosomal BED/BIM/FAM files...\n")
snpgdsBED2GDS(
  bed.fn = bed.fn,
  bim.fn = bim.fn,
  fam.fn = fam.fn,
  out.gdsfn = out.gdsfn
)
cat("Conversion complete: genotype_autosomes.gds created successfully.\n")
