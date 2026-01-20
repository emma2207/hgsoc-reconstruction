#!/usr/bin/env nextflow

process BAMMATCHER {
    conda "${params.conda}/bam-matcher"
    publishDir "${params.outdir}/bam-matcher", mode: 'copy'
    errorStrategy 'ignore'
    
    input:
        path(bam)
    
    output:
        path("${params.dataset}/*_output_report.txt")
    
    script:
        """
        set -euo pipefail

        output_location="${params.dataset}"
        mkdir -p \$output_location

        count=${#bam[@]}

        # Run bam-matcher over all pairs of bams (without repetition)
        for ((i = 0; i < count; i++))
        do
            for ((j = i; j < count; j++))
            do
                bam_1="${bams[i]}"
                bam_2="${bams[j]}"

                # Define output filenames based on sample names
                sample_1=\$(basename "\${bam_1}" .bam)
                sample_2=\$(basename "\${bam_2}" .bam)
                
                output_file="\${output_location}/\${sample_1}_vs_\${sample_2}_output_report.txt"

                # Run BAMMATCHER
                python ${params.BAMMATCHER}/bam-matcher.py \\
                    -B1 \${bam_1} \\
                    -B2 \${bam_2} \\
                    -o \${output_file}
            done
        done
        """
}