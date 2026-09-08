#!/usr/bin/env nextflow

process PSEUDOBULK_BAM_TO_FASTQ {
    conda "${params.conda}/ntsm"
    publishDir "${params.outdir}/pseudobulk_fastq", mode: 'copy'

    input:
    tuple val(sample_name), path(bam)

    output:
    tuple val(sample_name), path("${sample_name}.R1.fastq.gz"), path("${sample_name}.R2.fastq.gz")

    script:
    """
    set -euo pipefail

    if [ "${params.read_type}" = "single" ]
    then
        # For BAM-derived pseudobulks, reads are often unpaired/alignment-only.
        # Keep all reads in R1 for single-end mode and emit an empty R2 placeholder.
        samtools fastq -n \
            -0 ${sample_name}.R1.fastq \
            -1 /dev/null \
            -2 /dev/null \
            -s /dev/null \
            ${bam}
        : > ${sample_name}.R2.fastq
    else
        samtools fastq -n \
            -1 ${sample_name}.R1.fastq \
            -2 ${sample_name}.R2.fastq \
            -0 ${sample_name}.R0.fastq \
            -s ${sample_name}.RS.fastq \
            ${bam}

        # If BAM lacks read1/read2 flags, recover reads from unpaired outputs.
        if [ ! -s ${sample_name}.R1.fastq ] && [ ! -s ${sample_name}.R2.fastq ] && { [ -s ${sample_name}.R0.fastq ] || [ -s ${sample_name}.RS.fastq ]; }
        then
            cat ${sample_name}.R0.fastq ${sample_name}.RS.fastq > ${sample_name}.R1.fastq
            : > ${sample_name}.R2.fastq
        fi

        rm -f ${sample_name}.R0.fastq ${sample_name}.RS.fastq
    fi

    gzip -f ${sample_name}.R1.fastq
    gzip -f ${sample_name}.R2.fastq

    if [ "${params.read_type}" = "single" ]
    then
        if [ ! -s ${sample_name}.R1.fastq.gz ]
        then
            echo "ERROR: samtools fastq did not produce reads for ${sample_name} in single-end mode" >&2
            exit 1
        fi
    else
        if [ ! -s ${sample_name}.R1.fastq.gz ] && [ ! -s ${sample_name}.R2.fastq.gz ]
        then
            echo "ERROR: samtools fastq did not produce paired or recoverable reads for ${sample_name}" >&2
            exit 1
        fi
    fi
    """
}
