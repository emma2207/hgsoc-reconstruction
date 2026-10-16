#!/usr/bin/env nextflow

include { buildReadDepthContext } from './helpers/read_depth_utils'

process VIREO_MATCH {
    conda "${params.conda}/sample-matching"
    publishDir "${params.outdir}/2a_vireo", mode: 'copy'
    errorStrategy 'ignore'
    
    input:
        path(vcf_files)
        val(modalities)
    
    output:
        path("${params.dataset}/**/similarity_matrix.csv")
        path("${params.dataset}/**/matched_samples.csv")

    script:
        def readDepthContext = buildReadDepthContext(
            params.read_depth,
            params.dataset,
            params.pseudobulk,
            modalities,
            params.ncells
        )

        def outputPath = "${readDepthContext.experimentPath}/read_depth_${params.pseudobulk ? readDepthContext.pseudobulkReadDepth : readDepthContext.modalityReadDepthTag}"

        """
        set -euo pipefail

        if [ ${params.pseudobulk} == true ]
        then 
            output_location="${outputPath}"
            mkdir -p \$output_location

            python ${params.projectDir}/modules/vireo.py \
                -v1 ${vcf_files} \
                -v2 ${vcf_files} \
                -o \${output_location}
        else
            output_location="${outputPath}"
            mkdir -p \$output_location

            python ${params.projectDir}/modules/vireo.py \
                -v1 ${vcf_files[0]} \
                -v2 ${vcf_files[1]} \
                -o \${output_location}
        fi
        """
}