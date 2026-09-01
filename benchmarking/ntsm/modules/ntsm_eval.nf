#!/usr/bin/env nextflow

process NTSM_EVAL {
    conda "${params.conda}/ntsm"
    publishDir "${params.outdir}/ntsm_eval", mode: 'copy'

    input:
    path(count_files)

    output:
    path('ntsm_pairwise.tsv', optional: true)

    script:
    """
    set -euo pipefail

    ntsmEval ${count_files} > ntsm_pairwise.tsv
    """
}
