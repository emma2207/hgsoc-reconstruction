#!/usr/bin/env nextflow

process NGSCHECKMATE {
    tag "ngscheckmate"
    conda "${params.conda}/ngscheckmate"
    publishDir "${params.outdir}/ngscheckmate", mode: 'copy'
    errorStrategy 'ignore'
    
    input:
        path(vcf_path)
    
    output:
        tuple path("${params.dataset}_all.txt"), 
            path("${params.dataset}_matched.txt"), 
            path("${params.dataset}_output_corr_matrix.txt"), emit: results
    
    script:
        """
        set -euo pipefail

        echo \$(ls $vcf_path)
        
        # Run NGSCheckMate in VCF mode
        python ${params.NGS_CHECKMATE}/ncm.py \\
            -V \\
            -d ${vcf_path} \\
            -bed ${params.NGS_CHECKMATE}/SNP/SNP_GRCh38_hg38_wChr.bed \\
            -O . \\
            -N ${params.dataset}
        """
}