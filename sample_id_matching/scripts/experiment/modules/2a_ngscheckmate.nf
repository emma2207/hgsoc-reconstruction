#!/usr/bin/env nextflow

process NGSCHECKMATE {
    tag "ngscheckmate"
    conda "${params.conda}/ngscheckmate"
    publishDir "${params.outdir}/ngscheckmate", mode: 'copy'
    errorStrategy 'ignore'
    
    input:
        path(vcf_path)
    
    output:
        path("${params.dataset}/output_all.txt")
        path("${params.dataset}/output_matched.txt")
        path("${params.dataset}/output_corr_matrix.txt")
    
    script:
        """
        set -euo pipefail

        output_location="${params.dataset}"
        mkdir -p \$output_location

        echo \$(ls $vcf_path)
        
        # Run NGSCheckMate in VCF mode
        python ${params.NGS_CHECKMATE}/ncm.py \\
            -V \\
            -d ${vcf_path} \\
            -bed ${params.NGS_CHECKMATE}/SNP/SNP_GRCh38_hg38_wChr.bed \\
            -O \${output_location} 
        """
}