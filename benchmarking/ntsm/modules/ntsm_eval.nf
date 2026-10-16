#!/usr/bin/env nextflow

process NTSM_EVAL {
    conda "${params.conda}/ntsm"
    publishDir "${params.outdir}/ntsm_eval", mode: 'copy'

    input:
    path(count_files)
    val(modalities)

    output:
    path('ntsm_pairwise_*.tsv', optional: true)

    script:
    def modalitiesStr = (modalities instanceof Collection) ? modalities.join('_') : modalities.toString()
    def outputFile = params.pseudobulk
        ? "ntsm_pairwise_pseudobulk_${params.dataset}_ncells_${params.ncells}.tsv"
        : "ntsm_pairwise_${modalitiesStr}_${params.dataset}.tsv"
    """
    set -euo pipefail

    ntsmEval ${count_files} > ${outputFile}
    """
}
