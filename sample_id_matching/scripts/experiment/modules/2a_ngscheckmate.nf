#!/usr/bin/env nextflow

process NGSCHECKMATE {
    conda "${params.conda}/ngscheckmate"
    publishDir "${params.outdir}/2a_ngscheckmate", mode: 'copy'
    errorStrategy 'ignore'
    
    input:
        path(vcf)
    
    output:
        path("${params.dataset}/*/*/*/output_all.txt")
        path("${params.dataset}/*/*/*/output_matched.txt")
        path("${params.dataset}/*/*/*/output_output_corr_matrix.txt")

    script:
        """
        set -euo pipefail

        if [ ${params.pseudobulk} == true ]
        then 
            output_location="${params.dataset}/pseudobulk/ncells_${params.ncells}/read_depth_${params.read_depth}"
        else
            output_location="${params.dataset}/real_data/ncells_null/read_depth_${params.read_depth}" 
        fi
        mkdir -p \$output_location

        # Unzip and list
        for file in ${vcf}
        do  
            if ! [[ \$file == *"2507"* && \$file == *"bulk_diss_polyA"* ]]
            then
                base_name=\$(basename "\$file")
                base_name=\${base_name%.gz}

                gunzip -c \${file} > \${output_location}/\${base_name}
                
                echo "\${output_location}/\${base_name}" >> \${output_location}/vcf_list.txt
            fi
        done
        
        # Run NGSCheckMate in VCF mode
        python ${params.NGS_CHECKMATE}/ncm.py \\
            -V \\
            -l \${output_location}/vcf_list.txt \\
            -bed ${params.NGS_CHECKMATE}/SNP/SNP_GRCh38_hg38_wChr.bed \\
            -O \${output_location}
        """
}