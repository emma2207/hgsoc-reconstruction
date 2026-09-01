#!/usr/bin/env nextflow

process NTSM_COUNT {
    conda "${params.conda}/ntsm"
    publishDir "${params.outdir}/ntsm_counts", mode: 'copy'

    input:
    tuple val(sample_name), path(r1_fastqs), path(r2_fastqs)

    output:
    path("counts_${sample_name}.txt")

    script:
    """
    set -euo pipefail

    output_location="ntsm_counts/${params.dataset}/${params.datatype}"
    mkdir -p \${output_location}

    cat ${r1_fastqs} > \${output_location}/${sample_name}.R1.fastq.gz
    cat ${r2_fastqs} > \${output_location}/${sample_name}.R2.fastq.gz

    ntsmCount -t ${task.cpus} -s ${params.sites} \
        \${output_location}/${sample_name}.R1.fastq.gz \
        \${output_location}/${sample_name}.R2.fastq.gz > counts_${sample_name}.txt
    """
}
