## 🛠 GWAS of HMM-derived brain dynamics, with downstream analyses

| Script | Description |
|--------|-------------|
| **biobank_pipeline_hmm.m** | Matlab pipeline for HMM-based brain state analysis in UK Biobank. |
| **make_grm.sh** | Creates a genetic relationship matrix (GRM) for heritability and GWAS analysis. |
| **fastgwa.sh** | Runs fastGWA for rapid GWAS association testing. |
| **greml-heritability.sh** | Estimates SNP-based heritability using GCTA-GREML. |
| **ldsc-heritability.sh** | Estimates heritability using LD Score Regression. |
| **s-ldsc.sh** | Runs stratified LD Score Regression for partitioned heritability. |
| **clumping.sh** | Performs SNP clumping for GWAS summary statistics. |
| **extract_snp.sh** | Extracts specific SNPs from summary statistics for downstream analysis. |
| **munge.sh** | Prepares GWAS summary statistics for genetic correlation. |
| **genetic-corr.sh** | Computes genetic correlations between traits. |
| **magma.sh** | Runs MAGMA for positional gene mapping. |
| **h-magma.sh** | Runs H-MAGMA for chromatin-informed gene mapping. |
| **combo.sh** | Runs gene mapping considering masking effects. |
| **smr-query.R** | xQTL queries for gene mapping. |
| **smr.sh** | Runs SMR (Summary-based Mendelian Randomisation) for gene mapping. |
| **gwas_catalog | Links reported GWAS traits to the respective gene. |
| **polyfun1.sh** | First step in functionally informed fine-mapping. |
| **polyfun2.sh** | Second step in functionally informed fine-mapping. |
| **polyfun3.sh** | Third step in functionally informed fine-mapping. |
| **run_coloc.R** | Performs coloc.abf and coloc.susie analysis between GWAS traits. |
| **run_hyprcoloc.R** | Runs HyPrColoc for multi-trait colocalisation analysis. |
| **run_mr.R** | Performs Mendelian Randomisation analyses. |
| **run_mrlap.R** | Runs MRlap for Mendelian Randomisation with LD adjustment. |
