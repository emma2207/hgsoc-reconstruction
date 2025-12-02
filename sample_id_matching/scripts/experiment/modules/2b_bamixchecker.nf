#!/usr/bin/env nextflow

process BAMIXCHECKER {
    conda "${params.conda}/bamixchecker"
    publishDir "${params.outdir}/bamixchecker", mode: 'copy'
    errorStrategy 'ignore'
    cache false
    
    input:
        path(bam)
    
    output:
        path("${params.dataset}/BAMixChecker/BAMixChecker_Report.html")
        path("${params.dataset}/BAMixChecker/BAMixChecker_Heatmap.pdf")
        path("${params.dataset}/BAMixChecker/Total_result.txt")
        path("${params.dataset}/BAMixChecker/Matched_samples.txt"), optional: true
        path("${params.dataset}/BAMixChecker/Mismatched_samples.txt"), optional: true
    script:
        """
        set -euo pipefail

        # Create config file for BAMixChecker
        echo "GATK=${params.conda}/bamixchecker/bin/gatk" > BAMixChecker.config
        echo "BEDTOOLS=${params.conda}/bamixchecker/bin/bedtools" >> BAMixChecker.config

        # Fix the read groups
        for file in ${bam}
        do  
            base_name=\$(basename "\$file")
            base_name=\${base_name%.bam}

            ${params.conda}/bamixchecker/bin/picard AddOrReplaceReadGroups \\
                I=\$file \\
                O=\${base_name}_rg.bam \\
                SORT_ORDER=coordinate \\
                RGID=4 \\
                RGLB=lib1 \\
                RGPL=illumina \\
                RGPU=unit1 \\
                RGSM=pseudobulk_1 \\
                CREATE_INDEX=True

            echo "\${base_name}_rg.bam" >> bam_list.txt
        done

        # Run BAMixChecker
        python ${params.BAMIXCHECKER}/BAMixChecker.py \\
            -l bam_list.txt \\
            -r ${params.refGenome}/fasta/genome.fa \\
            -o ${params.dataset} \\
            -p 4
        """
}