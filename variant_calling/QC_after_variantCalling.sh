#!/bin/bash
# This script will run several QCs and
# generate some plots for sanity check.
# Gokberk Alagoz - 24.04.2025

# Set paths
inDir="/data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/results/variant_calling/filtered_variants"
outDir="/data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/results/variant_calling/post-variantCalling_QC"
logDir="${outDir}/logs"

# Create output dirs
mkdir -p "$outDir" "$logDir"

# Samples to process
samples=(
  "TAI3468"
  "TAI48406"
  "TAI48807"
  "TAI48879"
  "TAI661"
)

# Submit one job per sample
for SAMPLE in "${samples[@]}"; do
  INPUT_VCF="${inDir}/${SAMPLE}.filtered.snps.vcf.gz"

  # Create a job script for each sample
  JOB_SCRIPT="${logDir}/job_${SAMPLE}.sh"

    cat <<EOF > "$JOB_SCRIPT"
#!/bin/bash
#$ -N variantQC_${SAMPLE_ID}
#$ -cwd
#$ -q multi15.q
#$ -S /bin/bash
#$ -o ${logDir}/variantQC_${SAMPLE}.out     # STDOUT log
#$ -e ${logDir}/variantQC_${SAMPLE}.err     # STDERR log

module load bcftools/1.18 

echo "Processing sample ${SAMPLE}"

# Generate sumstats for per-sample
bcftools stats ${INPUT_VCF} > ${outDir}/${SAMPLE}.summary_stats.txt

EOF

# Submit the job to the grid
    qsub "$JOB_SCRIPT"

    echo "Submitted job for $SAMPLE"
done

echo "Submitted all jobs! Go have a coffee..."
