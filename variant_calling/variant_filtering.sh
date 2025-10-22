#/bin/bash
# This script will filter variants based on
# per-sample median sequencing depth values, 
# which are calculated using mosdepth, as well
# as hard-filter thresholds recommended by GATK
# best-practices guideline.
#
# Use bcftools to filter SNPs and indels separately.
# First, calculate min and max depth per-sample 
# based on median depths from mosdepth and min./max. depth
# calculations as described by Kuderna et al. (2023).
# min. depth = median_depth / 3, max. depth = median_depth x 2
#
# Then, filter SNPs and indels based on these custom 
# thresholds and hard thresholds recommended by GATK.
# for hard-filter thresholds, see the supplementary 
# material of Kuderna et al. (2023)
#
# Gokberk Alagoz - 22.04.2025

# Paths

inDir="/data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/results/variant_calling/genotypes"
outDir="/data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/results/variant_calling/filtered_variants"

mkdir -p "${OUT_DIR}"

# Sample info: SAMPLE MIN_DEPTH MAX_DEPTH INPUT_VCF
# You should hardcode the needed info here.

samples=(
  "TAI3468 11 66 ${inDir}/TAI3468_MKDN240009949-1A_E250024699_E250025014_L1.genotyped.vcf.gz"
  "TAI48406 4 24 ${inDir}/TAI48406_MKDN240009945-1A_E200024585_L1.genotyped.vcf.gz"
  "TAI48807 3 18 ${inDir}/TAI48807_MKDN240009946-1A_E200024585_L1.genotyped.vcf.gz"
  "TAI48879 6 38 ${inDir}/TAI48879_MKDN240009944-1A_E200024585_L1.genotyped.vcf.gz"
  "TAI661 11 68 ${inDir}/TAI661_MKDN240009947-1A_E200024585_E250025003_L1.genotyped.vcf.gz"
)

module load bcftools/1.18

for sample_line in "${samples[@]}"; do
  read -r SAMPLE MIN_DEPTH MAX_DEPTH INPUT_VCF <<< "${sample_line}"

  qsub -N filter_${SAMPLE} << EOF
#!/bin/bash
module load bcftools/1.18

echo "Starting bcftools filtering for sample: ${SAMPLE}"

# Filter SNPs
bcftools filter -e "TYPE!='snp' | (GT='het' & FMT/AD[*:*] < 3) | AC > 2 | FMT/DP <= ${MIN_DEPTH} | FMT/DP >= ${MAX_DEPTH} | QD < 2 | FS > 60 | MQ < 40 | SOR > 3 | ReadPosRankSum < -8.0 | MQRankSum < -12.5" -O z -o ${outDir}/${SAMPLE}.filtered.snps.vcf.gz ${INPUT_VCF}

# Filter Indels
bcftools filter -e "TYPE!='indel' | (GT='het' & FMT/AD[*:*] < 3) | FMT/DP <= ${MIN_DEPTH} | FMT/DP >= ${MAX_DEPTH} | QD < 2 | FS > 200 | MQ < 40 | SOR > 3 | ReadPosRankSum < -20.0" -O z -o ${outDir}/${SAMPLE}.filtered.indels.vcf.gz ${INPUT_VCF}

echo "Completed bcftools filtering for sample: ${SAMPLE}"
EOF

done
