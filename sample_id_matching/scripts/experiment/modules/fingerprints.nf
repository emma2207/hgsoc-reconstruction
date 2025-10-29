#!/usr/bin/env nextflow

process CROSSCHECK_FINGERPRINTS {
    tag "fingerprint_check"
    publishDir "${params.outdir}/fingerprints", mode: 'copy'
    
    input:
    tuple val(sample_id), path(vcf_file)
    
    output:
    path("fingerprint_matrix.txt"), emit: results
    
    script:
    """
    # Run Picard CrosscheckFingerprints
    picard CrosscheckFingerprints \
        I=${vcf_file} \
        O=fingerprint_matrix.txt \
        MATRIX_OUTPUT=true
    """

    stub:
    """
    touch fingerprint_matrix.txt
    """
}