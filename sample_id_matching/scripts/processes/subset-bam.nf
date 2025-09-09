#!/usr/bin/env nextflow

process SPLIT_READS {
    
    label 'split'

    input:
    tuple path(bam)

    output:
    path "split/*"

    stub:
    """
    echo "Splitting reads for sample ${params.sample} in dataset ${params.dataset}"
    """

    script:
    """
    mkdir -p "${params.outputDir}/subset-bam/${params.dataset}/${params.sample}"
    output_location="${params.outputDir}/subset-bam/${params.dataset}/${params.sample}"


    subset-bam -b ${bam} \
        -c ${params.petaLibrary}/${params.dataset}/${params.sample}/*.tsv \
        -o "${output_location}/subset.bam" \
        --cores 1
    """
}

barcodes_csv=$5