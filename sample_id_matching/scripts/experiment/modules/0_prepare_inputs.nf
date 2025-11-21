#!/usr/bin/env nextflow

process PREPARE_INPUTS {
    tag "prepare inputs"
    conda "${params.conda}/sample-matching"
    publishDir "${params.outdir}/prepared_bams", mode: 'copy'
    
    input:
        tuple val(sample_id), val(nr), path(bam_files)
    
    output:
        path("*.bam")

    script:
        """
        # Rename BAM files using the filename and the folder name
        # Original paths are "sample_id/pseudobulk_nr.bam"
        # The new names should be "sample_id_pseudobulk_nr.bam"
        for bam_file in ${bam_files}
        do
            new_name="${sample_id}_${nr}.bam"
            mv "\$bam_file" "\$new_name"

            echo "Renamed \$bam_file to \$new_name"
        done
        """
}