#!/usr/bin/env nextflow

process BAMIXCHECKER {
    tag "${sample_id}"
    conda "${params.conda}/bamixchecker"
    publishDir "${params.outdir}/bamixchecker", mode: 'copy'
    // errorStrategy 'ignore'
    
    input:
        path(bam)
    
    output:
        path("BAMixChecker_report.html")
        path("BAMixChecker_heatmap.pdf")
        path("Total_result.txt")
        path("Matched_samples.txt"), optional: true
        path("Mismatched_samples.txt"), optional: true

    script:
        """
        set -euo pipefail

        # Create config file for BAMixChecker
        echo "GATK=${params.conda}/bamixchecker/bin/gatk" > BAMixChecker.config
        echo "BEDTOOLS=${params.conda}/bamixchecker/bin/bedtools" >> BAMixChecker.config

        # Create input file list
        for file in ${bam}
        do
            echo "\$file" >> bam_list.txt
        done

        # Run BAMixChecker
        python ${params.BAMIXCHECKER}/BAMixChecker.py \\
          -l bam_list.txt \\
          -r ${params.refGenome}/fasta/genome.fa \\
          -o . \\
          -p 4
        """
}