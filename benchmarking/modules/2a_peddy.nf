#!/usr/bin/env nextflow

include { buildReadDepthContext } from './helpers/read_depth_utils'

process PEDDY {
    conda "${params.conda}/peddy"
    publishDir "${params.outdir}/2a_peddy", mode: 'copy'
    errorStrategy 'ignore'

    input:
        path(all_variants_vcf)
        path(all_variants_vcf_index)
        val(modalities)

    output:
        path("${params.dataset}/**/peddy*")

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

        if [ ${params.pseudobulk} == true ]
        then
            output_location="${params.dataset}/pseudobulk/ncells_${params.ncells}/read_depth_${pseudobulkReadDepth}"
        else
            output_location="${params.dataset}/${realDataNcellsPath}/read_depth_${modalityReadDepthTag}"
        fi
        mkdir -p \$output_location

        # Create a PED file for peddy analysis from VCF sample IDs.
        peddy_ped_file="\${output_location}/peddy.ped"
        bcftools query -l ${all_variants_vcf} | while IFS= read -r sample_id
        do
            [ -n "\$sample_id" ] || continue
            echo -e "\${sample_id}\t\${sample_id}\t0\t0\t0\t-9" >> "\$peddy_ped_file"
        done

        python -m peddy -p 4 \\
            --sites hg38 \\
            ${all_variants_vcf} \\
            \${peddy_ped_file} 
        """
}
