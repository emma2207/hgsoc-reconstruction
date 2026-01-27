#!/usr/bin/env nextflow

process BAMMATCHER {
    conda "${params.conda}/bam-matcher"
    publishDir "${params.outdir}/bam-matcher", mode: 'copy'
    errorStrategy 'ignore'
    
    input:
        tuple path(bam), path(bai)
    
    output:
        path("${params.dataset}/*_output_report.txt")
    
    script:
        """
        set -euo pipefail

        output_location="${params.dataset}"
        mkdir -p \$output_location

        # Convert paired bam/bai inputs to arrays
        bams=($bam)
        count=0
        for file in \$bams; do
            count=\$((count + 1))
        done

        # Process BAM files
        processed_bams=()
        for ((i = 0; i < count; i++))
        do
            original_bam="\${bams[i]}"
            sample_name=\$(basename "\${original_bam}" .bam)
            processed_bam="\${sample_name}_with_RG.bam"
            
            # Add read groups using samtools
            samtools addreplacerg \\
                -r ID:\${sample_name} \\
                -r LB:lib1 \\
                -r PL:ILLUMINA \\
                -r PU:unit1 \\
                -r SM:\${sample_name} \\
                -o temp_\${processed_bam} \\
                \${original_bam}
            
            processed_bams+=(\${processed_bam})
        done

        # Run bam-matcher over all pairs of processed bams (without repetition)
        for ((i = 0; i < count; i++))
        do
            for ((j = i; j < count; j++))
            do
                bam_1="\${processed_bams[i]}"
                bam_2="\${processed_bams[j]}"

                # Define output filenames based on original sample names  
                sample_1=\$(basename "\${bams[i]}" .bam)
                sample_2=\$(basename "\${bams[j]}" .bam)
                
                output_file="\${output_location}/\${sample_1}_vs_\${sample_2}_output_report.txt"

                # Clean cache directory
                if [ -d "\$output_location/cache" ]; then
                    rm -rf "\$output_location/cache"
                fi
                mkdir -p "\$output_location/cache"

                # Run BAMMATCHER
                echo "Running bam-matcher for \${sample_1} vs \${sample_2}"
                python ${params.BAMMATCHER}/bam-matcher.py \\
                    -B1 \${bam_1} \\
                    -B2 \${bam_2} \\
                    -o \${output_file} \\
                    --scratch-dir \$output_location/cache \\
                    2>&1 | tee bam_matcher_\${sample_1}_vs_\${sample_2}.log
            done
        done
        """
}