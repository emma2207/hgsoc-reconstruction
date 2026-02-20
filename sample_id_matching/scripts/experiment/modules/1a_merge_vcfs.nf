#!/usr/bin/env nextflow

process MERGE_AND_FILTER_VCFS {
    conda "${params.conda}/sample-matching"
    publishDir "${params.outdir}/1a_vcf", mode: 'copy'
    
    input:
        path individual_vcfs
        path index_files
        val modalities 

    
    output:
        path("${params.dataset}/*/*/*/*_modality_variants.vcf.gz"), emit: modality_vcfs
        path("${params.dataset}/*/*/*/*_individual_variants.vcf.gz"), emit: individual_vcfs
        path("${params.dataset}/*/*/*/all_variants.vcf.gz")
        path("${params.dataset}/*/*/*/filtered_variants_rd_*.vcf.gz")

    script:
        """
        if [ ${params.pseudobulk} == true ]
        then 
            output_location="${params.dataset}/pseudobulk/ncells_${params.ncells}/read_depth_${params.read_depth}"
        else
            output_location="${params.dataset}/real_data/ncells_null/read_depth_${params.read_depth}" 
        fi
        all_variants_output="\${output_location}/all_variants.vcf.gz"
        rd_filtered_variants="\${output_location}/filtered_variants_rd_${params.read_depth}.vcf.gz"
        mkdir -p \$output_location

        echo "Start merging individual VCFs"

        # Merge individual VCF files
        bcftools merge -Oz -o \${all_variants_output} ${individual_vcfs}

        echo "Finished merging VCFs!"

        # Index the merged VCF
        bcftools index \${all_variants_output}

        echo "Finished indexing..."

        # Filter variants further
        bcftools view -Oz -i 'DP>=${params.read_depth}' -o \${rd_filtered_variants} \${all_variants_output}
        bcftools index \${rd_filtered_variants}

        echo "Filtered variants by read depth"

        # List samples by modality
        for mod in ${modalities.join(' ')}
        do  
            echo "\$mod"
            for sample in `bcftools query -l \${rd_filtered_variants}`
            do
                if [[ "\$sample" == *"\$mod"* ]]
                then
                    echo "\$sample" >> "\${output_location}/\${mod}_files.txt"
                fi
            done
        done

        # Save variants separately for different modalities
        for mod in ${modalities.join(' ')}
        do
            bcftools view -S "\${output_location}/\${mod}_files.txt" \
                -Oz -o "\${output_location}/\${mod}_modality_variants.vcf.gz" \
                \${rd_filtered_variants}
            bcftools index \${output_location}/\${mod}_modality_variants.vcf.gz
        done

        echo "Saved variants per modality"

        # Filter individual VCF inputs by read depth
        for vcf_file in ${individual_vcfs}
        do
            sample_id=\$(basename "\$vcf_file" .vcf.gz)
            # sample_id="\${sample_name%????}"
            echo "Filtering \$sample_id"
            
            bcftools view \
                -c1 -Oz \
                -i "DP>=${params.read_depth}" \
                -o "\${output_location}/\${sample_id}_individual_variants.vcf.gz" \
                \$vcf_file

            bcftools index \${output_location}/\${sample_id}_individual_variants.vcf.gz
        done

        echo "Finished filtering individual samples by read depth"
        """
}