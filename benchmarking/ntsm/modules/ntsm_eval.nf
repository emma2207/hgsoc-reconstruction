#!/usr/bin/env nextflow

process NTSM_EVAL {
    conda "${params.conda}/ntsm"
    publishDir "${params.outdir}/ntsm_eval", mode: 'copy'

    input:
    path(count_files)
    val(modalities)

    output:
    path('ntsm_pairwise.tsv', optional: true)

    script:
    """
    set -euo pipefail

    if ${params.pseudobulk}
    then
        output_file="ntsm_pairwise_pseudobulk_${params.dataset}_ncells_${params.ncells}.tsv"
    else
        # Combine values of modalities into a single string for the output filename
        modalities_str = modalities.join('_')
        output_file="ntsm_pairwise_${modalities_str}_${params.dataset}.tsv"
    fi
    

    ntsmEval ${count_files} > ${output_file}
    """
}
