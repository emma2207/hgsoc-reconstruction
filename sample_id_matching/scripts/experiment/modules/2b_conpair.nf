#!/usr/bin/env nextflow

process CONPAIR {
    conda "${params.conda}/conpair"
    publishDir "${params.outdir}/conpair", mode: 'copy'
    errorStrategy 'ignore'
    
    input:
        path(bam)
    
    output:
        path("${params.dataset}/*.bam")
        path("${params.dataset}/*_concordance.txt")
    
    script:
        """
        set -euo pipefail

        export CONPAIR_DIR=/projects/${USER}/software/conpair
        export GATK_JAR=/projects/${USER}/software/anaconda/envs/conpair/opt/gatk-3.8/GenomeAnalysisTK.jar
        export PYTHONPATH=/projects/${USER}/software/conpair/modules

        output_location="${params.dataset}"
        mkdir -p \$output_location

        # Convert bam input to bash array and count files
        bams=($bam)
        count=0
        for file in $bam; do
            count=\$((count + 1))
        done

        # Run conpair over all pairs of bams (without repetition)
        for ((i = 0; i < count; i++))
        do  
            bam_1="\${bams[i]}"
            sample_1=\$(basename "\${bam_1}" .bam)

            python ${params.CONPAIR}/scripts/run_gatk_pileup_for_sample.py \\
                    -B \${bam_1} \\
                    -R ${params.refGenome}/fasta/genome.fa \\
                    -O \${output_location}/\${sample_1}.bam

            for ((j = i; j < count; j++))
            do
                bam_2="\${bams[j]}"
                sample_2=\$(basename "\${bam_2}" .bam)

                output_file="\${output_location}/\${sample_1}_vs_\${sample_2}_concordance.txt"
                
                # Run CONPAIR
                python ${params.CONPAIR}/scripts/run_gatk_pileup_for_sample.py \\
                    -B \${bam_2} \\
                    -R ${params.refGenome}/fasta/genome.fa \\
                    -O \${output_location}/\${sample_2}.bam

                python ${params.CONPAIR}/scripts/verify_concordance.py \\
                    -T \${sample_1}.bam \\
                    -N \${sample_2}.bam \\
                    -O \${output_file}
            done

        done
        """
}