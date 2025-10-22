#!/bin/bash
# This script will count variants in VCF files, prior to and after filtering.
# Gokberk Alagoz - 15.04.2025

# Paths

inDir="/data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/results/variant_calling/genotypes"
outDir="/data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/results/variant_calling/genotypes"

# Find all GVCF files
find "$inDir" -type f -name "*_L1.genotyped.vcf.gz" | while read -r VCF_FILE; do

    # Extract sample parameters
    ANIMAL_ID=$(basename "$VCF_FILE" | cut -d '_' -f 1)

    # Create a job script for each sample
    JOB_SCRIPT="${outDir}/logs/job_variantCount_${ANIMAL_ID}.sh"

    cat <<EOF > "$JOB_SCRIPT"
#!/bin/bash
#$ -N genotyping_${ANIMAL_ID}
#$ -cwd
#$ -q multi15.q
#$ -S /bin/bash
#$ -o ${outDir}/logs/${ANIMAL_ID}_variantCount.out     # STDOUT log
#$ -e ${outDir}/logs/${ANIMAL_ID}_variantCount.err     # STDERR log

module load java/jre1.8.0_201 gatk/4.3.0.0

echo "Processing sample: ${ANIMAL_ID}"

gatk CountVariants -V ${VCF_FILE}

echo "Finished processing ${ANIMAL_ID}"
EOF

    # Submit the job to the grid
    qsub "$JOB_SCRIPT"

    echo "Submitted job for $ANIMAL_ID"
done

echo "All jobs submitted!"
