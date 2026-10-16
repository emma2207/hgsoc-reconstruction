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
            modalities,
            params.ncells
        )

        def outputPath = "${readDepthContext.experimentPath}/read_depth_${params.pseudobulk ? readDepthContext.pseudobulkReadDepth : readDepthContext.modalityReadDepthTag}"

        """
        set -euo pipefail

        output_location="${outputPath}"
        mkdir -p \$output_location

        input_vcf=\$(realpath "${all_variants_vcf}")
        input_index=\$(realpath "${all_variants_vcf_index}")
        ln -sf "\$input_vcf" peddy.vcf.gz
        ln -sf "\$input_index" peddy.vcf.gz.csi
        peddy_vcf="\$PWD/peddy.vcf.gz"
        output_location="\$PWD/\$output_location"

        # Create a PED file for peddy analysis from VCF sample IDs.
        peddy_ped_file="\${output_location}/peddy.ped"
        bcftools query -l ${all_variants_vcf} | while IFS= read -r sample_id
        do
            [ -n "\$sample_id" ] || continue
            echo -e "\${sample_id}\t\${sample_id}\t0\t0\t0\t-9" >> "\$peddy_ped_file"
        done

        cd "\$output_location"
        python -m peddy -p 4 \\
            --sites hg38 \\
            "\$peddy_vcf" \\
            \${peddy_ped_file} 
        """
}
