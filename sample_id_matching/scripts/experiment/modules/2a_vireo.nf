#!/usr/bin/env nextflow

process VIREO_MATCH {
    conda "${params.conda}/sample-matching"
    publishDir "${params.outdir}/2a_vireo", mode: 'copy'
    errorStrategy 'ignore'
    
    input:
        path(vcf_files)
    
    output:
        path("${params.dataset}/*/similarity_matrix.csv")
        path("${params.dataset}/*/matched_samples.csv")
    
    script:
        """
        if [ ${params.pseudobulk} == true ]
        then 
            output_location="${params.dataset}/pseudobulk/ncells_${params.ncells}/read_depth_${params.read_depth}"
            mkdir -p \$output_location
            python ${params.projectDir}/modules/vireo.py \\
                -v1 ${vcf_files} \\
                -v2 ${vcf_files} \\
                -o ${params.dataset}/pseudobulk
        else
            output_location="${params.dataset}/real_data/read_depth_${params.read_depth}" 
            mkdir -p \$output_location

            python ${params.projectDir}/modules/vireo.py \\
                -v1 ${vcf_files[0]} \\
                -v2 ${vcf_files[1]} \\
                -o ${params.dataset}/real_data
        fi        
        """
}