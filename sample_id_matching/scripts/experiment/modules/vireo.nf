#!/usr/bin/env nextflow

process VIREO_MATCH {
    tag "vireo_matching"
    publishDir "${params.outdir}/vireo", mode: 'copy'
    
    input:
    tuple val(sample_id), path(vcf_file)
    
    output:
    path("vireo_results.txt"), emit: results
    
    script:
    """
    # Use the match_sample_id function from match_sample_ids.py
    python3 ${baseDir}/scripts/match_sample_ids.py \
        --vcf ${vcf_file} \
        --output vireo_results.txt
    """

    stub:
    """
    touch vireo_results.txt
    """
}