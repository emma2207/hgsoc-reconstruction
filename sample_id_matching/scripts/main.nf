#!/usr/bin/env nextflow

/*
 * main.nf — run simulated data pipeline
 */

nextflow.enable.dsl = 2

process QC_READS_WITH_FASTP {

    label '01_fastp'

    conda 'mismatch_project-install'

    input:
    tuple path(fastq)

    output:
    path "fastp/*"

    script:
    """
    bash scripts/run_fastp.sh ${fastq}
    """
}

process ALIGNMENT_WITH_STAR {

    label '02_star'

    conda 'mismatch_project-install'

    input:
    tuple path(fastq)

    output:
    path "star/*"

    script:
    """
    bash scripts/02_star_alignment.sh ${pool} ${fastq}
    """
}