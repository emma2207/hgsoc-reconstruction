#!/usr/bin/env nextflow

process NGSCHECKMATE {
    conda "${params.conda}/ngscheckmate"
    publishDir "${params.outdir}/ngscheckmate", mode: 'copy'
    errorStrategy 'ignore'
    cache false
    
    input:
        path(vcf)
    
    output:
        path("${params.dataset}/vcf_list.txt")
        path("${params.dataset}/output_all.txt")
        path("${params.dataset}/output_matched.txt")
        path("${params.dataset}/output_output_corr_matrix.txt")
    
    script:
        """
        set -euo pipefail

        output_location="${params.dataset}"
        mkdir -p \$output_location

        # Unzip and list
        for file in ${vcf}
        do  
            base_name=\$(basename "\$file")
            base_name=\${base_name%.gz}

            gunzip -c \${file} > \${output_location}/\${base_name}
            
            echo "\${output_location}/\${base_name}" >> \${output_location}/vcf_list.txt
        done
        
        # Run NGSCheckMate in VCF mode
        python ${params.NGS_CHECKMATE}/ncm.py \\
            -V \\
            -l \${output_location}/vcf_list.txt \\
            -bed ${params.NGS_CHECKMATE}/SNP/SNP_GRCh38_hg38_wChr.bed \\
            -O \${output_location}
        """
}