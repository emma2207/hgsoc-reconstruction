#!/usr/bin/env nextflow

process CREATE_PSEUDOBULK {
    
    label 'pseudobulk'

    input:
    tuple path(bam)

    output:
    path "pseudobulk/*"

    stub:
    """
    echo "Creating pseudobulk for sample ${params.sample} in dataset ${params.dataset}"
    """

    script:
    """
    mkdir -p "${params.outputDir}/pseudobulk/${params.dataset}/${params.sample}"
    output_location="${params.outputDir}/pseudobulk/${params.dataset}/${params.sample}"

    # Run pseudobulk creation
    samtools # Something something pseudobulk
    """
}