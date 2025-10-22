/*
    This script will follow the GATK pipeline for DNA variant calling.
    Gokberk Alagoz - 13.03.25
*/

/*
    Paths
*/

params.inDir = "/data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/results/mapping"
params.outDir = "/data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/results/variant_calling"
params.reGenome = "/data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/resources/genomes/chimpanzee/NHGRI_mPanTro3-v2.0_pri_ncbi/ncbi_dataset/data/GCF_028858775.2/GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.fna"
params.genomeIndex = "/data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/resources/genomes/chimpanzee/NHGRI_mPanTro3-v2.0_pri_ncbi/ncbi_dataset/data/GCF_028858775.2/panTro3-v2.0"

// Define the workflow
workflow {

    // Create the output directory
    process createOutDir {
        script:        
        """
        mkdir "${outDir}"
        """
    }

    // Step 1: Mark duplicates
    process markDuplicates {
        input:
        path bamFile from find("${inDir}", "*.sorted.bam").view()

        // Mark duplicates using GATK
        output:
        path "${outDir}/marked_${basename(bamFile)}.bam"

        script:
        """
        module load java/jdk1.8.0_201 gatk/4.3.0.0
        gatk MarkDuplicates -I ${bamFile} -O ${outDir}/marked_${basename(bamFile)}.bam -M ${outDir}/metrics_${basename(bamFile)}.txt
        """
    }

    // Step 2: Base quality score recalibration
    process baseRecalibration {
        input:
        path markedBamFile from find("${outDir}", "marked_*.bam")

        // Recalibrate base quality scores using GATK
        output:
        path "${outDir}/recalibrated_${basename(markedBamFile)}.bam"

        script:
        """
        module load java/jdk1.8.0_201 gatk/4.3.0.0
        gatk BaseRecalibrator -I ${markedBamFile} -knownSites ${refGenome} -O ${outDir}/recalibrated_${basename(markedBamFile)}.bam
        """
    }

    // Step 3: Variant calling
    process variantCalling {
        input:
        path recalibratedBamFile from find("${outDir}", "recalibrated_*.bam")

        // Call variants using GATK
        output:
        path "${outDir}/variants_${basename(recalibratedBamFile)}.vcf"

        script:
        """
        module load java/jdk1.8.0_201 gatk/4.3.0.0
        gatk HaplotypeCaller -I ${recalibratedBamFile} -O ${outDir}/variants_${basename(recalibratedBamFile)}.vcf -ERC GVCF
        """
    }
}
