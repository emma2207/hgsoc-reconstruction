#!/usr/bin/env nextflow

process NGSCHECKMATE {
    tag "ngscheckmate"
    publishDir "${params.outdir}/ngscheckmate", mode: 'copy'
    
    input:
    tuple val(sample_id), path(vcf_file)
    
    output:
    path("ngscheckmate_output.txt"), emit: results
    
    script:
    """
    # Run NGSCheckMate
    # Note: Adjust the command based on your NGSCheckMate installation and requirements
    ngscheckmate \
        --vcf ${vcf_file} \
        --outfile ngscheckmate_output.txt
    """

    stub:
    """
    touch ngscheckmate_output.txt
    """
}