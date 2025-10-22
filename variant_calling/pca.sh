#!/bin/bash
# This script will prepare filtered VCF files
# for PCA and then will perform PCA.
# Gokberk Alagoz - 30.04.2025

# Paths
inDir="/data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/results/variant_calling/filtered_variants"
outDir="/data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/results/variant_calling/PCA"

# Merge per-sample VCF files

bcftools merge TAI3468.filtered.snps.vcf.gz TAI48406.filtered.snps.vcf.gz TAI48807.filtered.snps.vcf.gz TAI48879.filtered.snps.vcf.gz TAI661.filtered.snps.vcf.gz -Oz -o merged.filtered.snps.vcf.gz
bcftools index merged.filtered.snps.vcf.gz

# Filtered VCFs contain multiallelic SNPs, which
# PLINK cannot work with. Remove these.

bcftools view -m2 -M2 -v snps merged.filtered.snps.vcf.gz -Oz -o merged.filtered.snps.biallelic.vcf.gz
bcftools index merged.filtered.snps.biallelic.vcf.gz

# Convert VCF files to PLINK files
# Make sure to add --allow-extra-chr option
# because PLINK does not recognize chimpanzee
# chromosome names.

plink2 --vcf merged.filtered.snps.biallelic.vcf.gz --make-bed --allow-extra-chr --out merged_filtered_snps

# Exctract principal components

plink2 --bfile merged_filtered_snps --pca --allow-extra-chr --out pca_results # errored out, requires N>50
