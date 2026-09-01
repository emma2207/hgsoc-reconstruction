#!/usr/bin/env nextflow

process PROCESS_METADATA {
    label 'process_metadata'
    conda "${params.conda}/ntsm"
    publishDir "${params.outdir}/metadata", mode: 'copy'

    input:
    path(metadata)
    val(fastq_dir)

    output:
    path('run_name_mapping.csv')

    script:
    """
    python ${params.projectDir}/process_metadata.py \
        --metadata \"$metadata\" \
        --fastq-dir \"$fastq_dir\" \
        --read-type ${params.read_type} \
        --output run_name_mapping.csv
    """
}
