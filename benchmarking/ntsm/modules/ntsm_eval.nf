#!/usr/bin/env nextflow

process NTSM_EVAL {
    conda "${params.conda}/ntsm"
    publishDir "${params.outdir}/ntsm_eval", mode: 'copy'

    input:
    path(count_files)

    output:
    path('ntsm_summary.txt')
    path('ntsm_pairwise.tsv', optional: true)

    script:
    def count_args = count_files.join(' ')
    def n_files = count_files.size()

    """
    set -euo pipefail

    echo "NTSM count files provided: ${n_files}" > ntsm_summary.txt

    ntsmEval ${count_files} > ntsm_pairwise.tsv
    """
}
