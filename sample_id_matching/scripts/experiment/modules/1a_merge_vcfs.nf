#!/usr/bin/env nextflow

process MERGE_AND_FILTER_VCFS {
    conda "${params.conda}/sample-matching"
    publishDir "${params.outdir}/1a_vcf", mode: 'copy'

    input:
        path individual_vcfs
        path index_files
        val modalities
        val is_pseudobulk

    output:
        path("${params.dataset}/*/*/*/*_modality_variants.vcf.gz"), emit: modality_vcfs
        path("${params.dataset}/*/*/*/*_individual_variants.vcf.gz"), emit: individual_vcfs
        path("${params.dataset}/*/*/*/all_variants.vcf.gz")
        path("${params.dataset}/*/*/*/filtered_variants_*_rd_*.vcf.gz")

    script:
            def resolvedReadDepth =
                (params.read_depth instanceof Map)
                    ? (params.read_depth[params.dataset] ?: params.read_depth.default ?: ["default": 0])
                    : ["default": params.read_depth]

            def resolvedDefault =
                (resolvedReadDepth instanceof Map)
                    ? (resolvedReadDepth.default ?: 0)
                    : resolvedReadDepth

            // For pseudobulk, force all modalities to use default read depth
            def datasetReadDepth =
                is_pseudobulk
                    ? ["default": resolvedDefault]
                    : resolvedReadDepth

            def readDepthPairs = datasetReadDepth.collect { k, v -> "${k}:${v}" }.join(',')

            def pseudobulkReadDepth =
                (datasetReadDepth instanceof Map)
                    ? (datasetReadDepth.default ?: 0)
                    : datasetReadDepth

        """
        if [ ${is_pseudobulk} == true ]
        then
            output_location="${params.dataset}/pseudobulk/ncells_${params.ncells}/read_depth_${pseudobulkReadDepth}"
        else
            output_location="${params.dataset}/real_data/ncells_null/read_depth_modality_specific"
        fi

        all_variants_output="\${output_location}/all_variants.vcf.gz"
        mkdir -p "\$output_location"

        READ_DEPTH_MAP='${readDepthPairs}'
        MODALITIES_STR='${modalities.join(' ')}'

        declare -A RD_BY_MOD
        IFS=',' read -ra RD_PAIRS <<< "\$READ_DEPTH_MAP"
        for pair in "\${RD_PAIRS[@]}"
        do
            key="\${pair%%:*}"
            value="\${pair#*:}"
            RD_BY_MOD["\$key"]="\$value"
        done
        DEFAULT_RD="\${RD_BY_MOD[default]:-0}"

        read -r -a MODS <<< "\$MODALITIES_STR"

        echo "Start merging individual VCFs"

        # Merge individual VCF files
        bcftools merge -Oz -o "\${all_variants_output}" ${individual_vcfs}
        bcftools index "\${all_variants_output}"

        echo "Finished merging and indexing all variants"

        # Build sample lists by modality from merged VCF
        for mod in "\${MODS[@]}"
        do
            : > "\${output_location}/\${mod}_files.txt"
            for sample in \$(bcftools query -l "\${all_variants_output}")
            do
                if [[ "\$sample" == *"\$mod"* ]]
                then
                    echo "\$sample" >> "\${output_location}/\${mod}_files.txt"
                fi
            done
        done

        # Save variants separately for each modality with modality-specific read depth
        for mod in "\${MODS[@]}"
        do
            mod_rd="\${RD_BY_MOD[\$mod]:-\$DEFAULT_RD}"
            mod_filtered_variants="\${output_location}/filtered_variants_\${mod}_rd_\${mod_rd}.vcf.gz"

            bcftools view -Oz -i "DP>=\${mod_rd}" -o "\${mod_filtered_variants}" "\${all_variants_output}"
            bcftools index "\${mod_filtered_variants}"

            bcftools view -S "\${output_location}/\${mod}_files.txt" \
                -Oz -o "\${output_location}/\${mod}_modality_variants.vcf.gz" \
                "\${mod_filtered_variants}"
            bcftools index "\${output_location}/\${mod}_modality_variants.vcf.gz"
        done

        echo "Saved modality VCFs with modality-specific read depth filtering"

        # Filter individual VCF inputs by modality-specific read depth
        for vcf_file in ${individual_vcfs}
        do
            sample_id=\$(basename "\$vcf_file" .vcf.gz)
            sample_mod="default"

            for mod in "\${MODS[@]}"
            do
                if [[ "\$sample_id" == *"\$mod"* ]]
                then
                    sample_mod="\$mod"
                    break
                fi
            done

            sample_rd="\${RD_BY_MOD[\$sample_mod]:-\$DEFAULT_RD}"

            echo "Filtering \$sample_id (modality=\$sample_mod, read_depth=\$sample_rd)"

            bcftools view \
                -c1 -Oz \
                -i "DP>=\${sample_rd}" \
                -o "\${output_location}/\${sample_id}_individual_variants.vcf.gz" \
                "\$vcf_file"

            bcftools index "\${output_location}/\${sample_id}_individual_variants.vcf.gz"
        done

        echo "Finished filtering individual samples by modality-specific read depth"
        """
}