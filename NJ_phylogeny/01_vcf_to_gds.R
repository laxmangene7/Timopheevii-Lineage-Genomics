# ==============================================================================
# 01_vcf_to_gds.R
#
# Converts the filtered T. timopheevii SNP VCF to SNPRelate GDS format
# (biallelic SNPs only), then prints the sample and SNP counts and checks that
# the outgroup accessions are present.
#
# Input : timopheevii_WGS_SNPs.vcf.gz
# Output: timopheevii_WGS_SNPs.gds
# Usage : Rscript 01_vcf_to_gds.R
# ==============================================================================

suppressPackageStartupMessages({
    library(gdsfmt)
    library(SNPRelate)
})

vcf_file <- "timopheevii_WGS_SNPs.vcf.gz"
gds_file <- "timopheevii_WGS_SNPs.gds"

# Wild emmer (T. dicoccoides) accessions used as the outgroup in step 2
outgroup_samples <- c("PI538450", "IG116176", "KU-1982", "TA873_Hakan")

if (!file.exists(vcf_file)) stop("VCF file not found: ", vcf_file)

message("Converting VCF to GDS (biallelic SNPs)...")
snpgdsVCF2GDS(vcf_file, gds_file, method = "biallelic")

gds <- snpgdsOpen(gds_file)
ids <- read.gdsn(index.gdsn(gds, "sample.id"))
n_snps <- length(read.gdsn(index.gdsn(gds, "snp.id")))
snpgdsClose(gds)

message("Samples: ", length(ids), " | SNPs: ", n_snps)

missing_og <- setdiff(outgroup_samples, ids)
if (length(missing_og) > 0)
    stop("Outgroup samples not found in the VCF: ", paste(missing_og, collapse = ", "))

message("Done: ", gds_file)
