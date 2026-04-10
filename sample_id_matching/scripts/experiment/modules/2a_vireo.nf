#!/usr/bin/env nextflow

process VIREO_MATCH {
    conda "${params.conda}/sample-matching"
    publishDir "${params.outdir}/2a_vireo", mode: 'copy'
    errorStrategy 'ignore'
    
    input:
        path(vcf_files)
    
    output:
        path("${params.dataset}/*/*/*/similarity_matrix.csv")
        path("${params.dataset}/*/*/*/matched_samples.csv")

    script:
        def resolvedReadDepth =
            (params.read_depth instanceof Map)
                ? (params.read_depth[params.dataset] ?: params.read_depth.default ?: ["default": 0])
                : ["default": params.read_depth]

        def pseudobulkReadDepth =
            (resolvedReadDepth instanceof Map)
                ? (resolvedReadDepth.default ?: 0)
                : resolvedReadDepth

        """
        set -euo pipefail

        if [ ${params.pseudobulk} == true ]
        then 
            output_location="${params.dataset}/pseudobulk/ncells_${params.ncells}/read_depth_${pseudobulkReadDepth}"
            mkdir -p \$output_location

            python ${params.projectDir}/modules/vireo.py \
                -v1 ${vcf_files} \
                -v2 ${vcf_files} \
                -o \${output_location}
        else
            output_location="${params.dataset}/real_data/ncells_null/read_depth_modality_specific"
            mkdir -p \$output_location

            python ${params.projectDir}/modules/vireo.py \
                -v1 ${vcf_files[0]} \
                -v2 ${vcf_files[1]} \
                -o \${output_location}
        fi
        """
}