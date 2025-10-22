#!/bin/bash
# This script will call haplotypes.
# Gokberk Alagoz - 07.04.2025

# Paths

inDir="/data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/results/mapping"
refFile="/data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/resources/genomes/chimpanzee/NHGRI_mPanTro3-v2.0_pri_ncbi/ncbi_dataset/data/GCF_028858775.2/GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.fna"
outDir="/data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/results/mapping/haplotype_calling"

# Make directories
mkdir -p "$outDir"
mkdir -p "${outDir}/logs"

# Find all read group-added BAM files
find "$inDir" -type f -name "*.addRG.marked.sorted.bam" | while read -r BAM_FILE; do
    # Extract sample paramaeters
    SAMPLE_ID=$(basename "$BAM_FILE" | sed 's/.addRG.marked.sorted.bam//')
    SAMPLE_RGID=$(basename "$BAM_FILE" | sed -r 's/.*1A_(.*?).addRG.marked.*/\1/')
    LIBRARY_ID=$(basename "$BAM_FILE" | cut -d '_' -f 2)
    ANIMAL_ID=$(basename "$BAM_FILE" | cut -d '_' -f 1)

    # Define output files
    BAM_sorted="${outDir}/${SAMPLE_ID}.addRG.marked.REsorted.bam"
    HC="${outDir}/${SAMPLE_ID}.raw.snps.indels.g.vcf.gz"

    # Create a job script for each sample
    JOB_SCRIPT="${outDir}/logs/job_${SAMPLE_ID}.sh"

    cat <<EOF > "$JOB_SCRIPT"
#!/bin/bash
#$ -N haplotype_calling_${SAMPLE_ID}
#$ -cwd
#$ -q multi15.q
#$ -S /bin/bash
#$ -o ${outDir}/logs/${SAMPLE_ID}.out     # STDOUT log
#$ -e ${outDir}/logs/${SAMPLE_ID}.err     # STDERR log

module load java/jre1.8.0_201 gatk/4.3.0.0 samtools/1.9

echo "Processing sample: $SAMPLE_ID"

# Sort and index final BAM files, just to 
# make sure everything is in order.
#samtools sort ${BAM_FILE} -o ${BAM_sorted}
#samtools index ${BAM_sorted}

gatk HaplotypeCaller -I ${BAM_sorted} \
                     -O ${HC} \
		     -R ${refFile} \
                     -ERC GVCF \

echo "Finished processing $SAMPLE_ID"
EOF

    # Submit the job to the grid
    qsub "$JOB_SCRIPT"

    echo "Submitted job for $SAMPLE_ID"
done

echo "All jobs submitted!"
