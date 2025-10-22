#!/bin/bash
# This script will generate callability masks that define
# which positions in the genome are reliably callable per
# sample/individual. This step is not directly related to
# variant quality, but about if the site is reliably
# sequenced and genotype-able.
# These masks are required to estimate heterozygosity, to
# define the callable genome size.
# Gokberk Alagoz - 23.04.2025

# Paths

inDir="/data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/results/variant_calling/genotypes"
outDir="/data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/results/variant_calling/callability_masks"

mkdir -p "${outDir}"
mkdir -p "${outDir}/logs"

# Generate callability masks

module load bcftools/1.18

# Sample information: sample MIN_HET_AD MIN_COV MAX_COV INPUT_VCF
samples=(
  "TAI3468 3 11 66 ${inDir}/TAI3468_MKDN240009949-1A_E250024699_E250025014_L1.genotyped.vcf.gz"
  "TAI48406 3 4 24 ${inDir}/TAI48406_MKDN240009945-1A_E200024585_L1.genotyped.vcf.gz"
  "TAI48807 3 3 18 ${inDir}/TAI48807_MKDN240009946-1A_E200024585_L1.genotyped.vcf.gz"
  "TAI48879 3 6 38 ${inDir}/TAI48879_MKDN240009944-1A_E200024585_L1.genotyped.vcf.gz"
  "TAI661 3 11 68 ${inDir}/TAI661_MKDN240009947-1A_E200024585_E250025003_L1.genotyped.vcf.gz"
)

for sample_line in "${samples[@]}"; do
  read -r SAMPLE MIN_HET_AD MIN_COV MAX_COV INPUT_VCF <<< "${sample_line}"

  qsub -N callmask_${SAMPLE} -cwd -V \
       -o "${outDir}/logs/${SAMPLE}.o" \
       -e "${outDir}/logs/${SAMPLE}.e" <<EOF
#!/bin/bash
module load bcftools/1.18

echo "Running callability mask filtering for sample: ${SAMPLE}"

bcftools filter -e "(GT='./.') | (GT='het' & FMT/AD[*:*] < ${MIN_HET_AD} ) | FMT/DP <= ${MIN_COV} | FMT/DP >= ${MAX_COV} | FMT/GQ <= 30" \\
  -O z -o ${outDir}/${SAMPLE}.callability_mask.vcf.gz \\
  ${INPUT_VCF}

echo "Finished callability mask filtering for sample: ${SAMPLE}"
EOF
done

