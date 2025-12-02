#!/usr/bin/env nextflow

process VIREO_MATCH {
    conda "${params.conda}/sample-matching"
    publishDir "${params.outdir}/vireo", mode: 'copy'
    errorStrategy 'ignore'
    cache false
    
    input:
        path(vcf_file)
    
    output:
        path("${params.dataset}/similarity_matrix.csv")
        path("${params.dataset}/matched_samples.csv")
    
    script:
        """
        output_location="${params.dataset}"
        mkdir -p \$output_location

        python ${params.projectDir}/modules/vireo.py \\
            -v ${vcf_file} \\
            -d ${params.dataset}
        """
}