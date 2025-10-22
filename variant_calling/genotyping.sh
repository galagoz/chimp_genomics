#!/bin/bash
# This script will genotype raw haplotype calls.
# Gokberk Alagoz - 10.04.2025

# Paths

inDir="/data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/results/variant_calling/haplotype_calling"
refFile="/data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/resources/genomes/chimpanzee/NHGRI_mPanTro3-v2.0_pri_ncbi/ncbi_dataset/data/GCF_028858775.2/GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.fna"
outDir="/data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/results/variant_calling/genotypes"

# Make directories
mkdir -p "$outDir"
mkdir -p "${outDir}/logs"

# Find all GVCF files
find "$inDir" -type f -name "*.raw.snps.indels.g.vcf.gz" | while read -r GVCF_FILE; do
    # Extract sample paramaeters
    SAMPLE_ID=$(basename "$GVCF_FILE" | sed 's/.raw.snps.indels.g.vcf.gz//')
    LIBRARY_ID=$(basename "$GVCF_FILE" | cut -d '_' -f 2)
    ANIMAL_ID=$(basename "$GVCF_FILE" | cut -d '_' -f 1)

    # Define output files
    GENOTYPE="${outDir}/${SAMPLE_ID}.genotyped.vcf.gz"

    # Create a job script for each sample
    JOB_SCRIPT="${outDir}/logs/job_${SAMPLE_ID}.sh"

    cat <<EOF > "$JOB_SCRIPT"
#!/bin/bash
#$ -N genotyping_${SAMPLE_ID}
#$ -cwd
#$ -q multi15.q
#$ -S /bin/bash
#$ -o ${outDir}/logs/${SAMPLE_ID}.out     # STDOUT log
#$ -e ${outDir}/logs/${SAMPLE_ID}.err     # STDERR log

module load java/jre1.8.0_201 gatk/4.3.0.0

echo "Processing sample: $SAMPLE_ID"

gatk GenotypeGVCFs --include-non-variant-sites \
                   -R ${refFile} \
                   -V ${GVCF_FILE} \
		   -O ${GENOTYPE}

echo "Finished processing $SAMPLE_ID"
EOF

    # Submit the job to the grid
    qsub "$JOB_SCRIPT"

    echo "Submitted job for $SAMPLE_ID"
done

echo "All jobs submitted!"
