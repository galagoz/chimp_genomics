#!/bin/bash
# This script will combine per-chimpanzee GVCF files
# produced by HaplotypeCaller into a single multi-sample GVCF.
# For details: https://gatk.broadinstitute.org/hc/en-us/articles/360037053272-CombineGVCFs
# Gokberk Alagoz - 10.04.2025

# Paths

inDir="/data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/results/variant_calling/haplotype_calling"
refFile="/data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/resources/genomes/chimpanzee/NHGRI_mPanTro3-v2.0_pri_ncbi/ncbi_dataset/data/GCF_028858775.2/GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.fna"
outDir="/data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/results/variant_calling/genotypes"

# Make directories
mkdir -p "$outDir"
mkdir -p "${outDir}/logs"

# Get list of all GVCF files to combine
GVCF_FILES=($(find "$inDir" -type f -name "*.raw.snps.indels.g.vcf.gz"))

# Build --variant arguments
VARIANT_ARGS=""
for GVCF in "${GVCF_FILES[@]}"; do
  VARIANT_ARGS+=" --variant $GVCF"
done

# Make a script to combine them all and in the darkness bind them.
JOB_SCRIPT="${outDir}/logs/job_combine_gvcfs.sh"

cat <<EOF > "$JOB_SCRIPT"
#!/bin/bash
#$ -N CombineGVCFs_and_genotype
#$ -cwd
#$ -q multi15.q
#$ -S /bin/bash
#$ -o ${outDir}/logs/combine_gvcfs_and_genotype.out   # STDOUT log
#$ -e ${outDir}/logs/combine_gvcfs_and_genotype.err   # STDERR log

module load java/jre1.8.0_201 gatk/4.3.0.0

gatk CombineGVCFs -R ${refFile} \
		  ${VARIANT_ARGS} \
		  -O "${outDir}/TAI_pilot_chimpanzees.vcf.gz"

EOF

qsub "$JOB_SCRIPT"
echo "Submitted the job!"
