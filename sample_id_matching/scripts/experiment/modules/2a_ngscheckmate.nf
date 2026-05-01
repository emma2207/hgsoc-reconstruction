#!/usr/bin/env nextflow

include { buildReadDepthContext } from './helpers/read_depth_utils'

process NGSCHECKMATE {
    conda "${params.conda}/ngscheckmate"
    publishDir "${params.outdir}/2a_ngscheckmate", mode: 'copy'
    errorStrategy 'ignore'
    
    input:
        path(vcf)
        val(modalities)
    
    output:
        path("${params.dataset}/**/output_all.txt")
        path("${params.dataset}/**/output_matched.txt")
        path("${params.dataset}/**/output_output_corr_matrix.txt")

    script:
        def readDepthContext = buildReadDepthContext(
            params.read_depth,
            params.dataset,
            params.pseudobulk,
            modalities
        )

        def modalityReadDepthTag = readDepthContext.modalityReadDepthTag
        def pseudobulkReadDepth = readDepthContext.pseudobulkReadDepth
        def realDataNcellsPath = readDepthContext.realDataNcellsPath

        """
        set -euo pipefail
z
        if [ ${params.pseudobulk} == true ]
        then 
            output_location="${params.dataset}/pseudobulk/ncells_${params.ncells}/read_depth_${pseudobulkReadDepth}"
        else
            output_location="${params.dataset}/${realDataNcellsPath}/read_depth_${modalityReadDepthTag}" 
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