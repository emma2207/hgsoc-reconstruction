#!/usr/bin/env nextflow

process CONPAIR {
    conda "${params.conda}/conpair"
    publishDir "${params.outdir}/conpair", mode: 'copy'
    errorStrategy 'ignore'
    
    input:
        path(bam)
    
    output:
        path("${params.dataset}/*_concordance.txt")
    
    script:
        """
        set -euo pipefail

        output_location="${params.dataset}"
        mkdir -p \$output_location

        count=${#bam[@]}

        # Run conpair over all pairs of bams (without repetition)
        for ((i = 0; i < count; i++))
        do
            for ((j = i; j < count; j++))
            do
                bam_1="${bams[i]}"
                bam_2="${bams[j]}"

                # Define output filenames based on sample names
                sample_1=\$(basename "\${bam_1}" .bam)
                sample_2=\$(basename "\${bam_2}" .bam)
                
                output_file="\${output_location}/\${sample_1}_vs_\${sample_2}_concordance.txt"

                # Run CONPAIR
                python ${params.CONPAIR}/scripts/verify_concordance.py \\
                    -T \${bam_1} \\
                    -N \${bam_2} \\
                    -O \${output_file}
            done

        done
        """
}