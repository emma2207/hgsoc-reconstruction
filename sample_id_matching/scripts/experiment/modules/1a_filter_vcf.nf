#!/usr/bin/env nextflow

process FILTER_VCF {
    conda "${params.conda}/sample-matching"
    publishDir "${params.outdir}/1a_vcf", mode: 'copy'
    
    input:
        path(vcf_file)
        val(modalities)

    
    output:
        path("${params.dataset}/*/*/*/*_individual_variants.vcf.gz"), emit: individual_vcfs
        path("${params.dataset}/*/*/*/*_modality_variants.vcf.gz")
        path("${params.dataset}/*/*/*/filtered_variants_rd_*.vcf.gz"), emit: modality_vcfs

    script:
        """
        if [ ${params.pseudobulk} == true ]
        then 
            output_location="${params.dataset}/pseudobulk/ncells_${params.ncells}/read_depth_${params.read_depth}"
        else
            output_location="${params.dataset}/real_data/ncells_null/read_depth_${params.read_depth}" 
        fi
        mkdir -p \$output_location

        # List samples by modality
        for sample in `bcftools query -l $vcf_file`
        do  
            for mod in ${modalities.join(' ')}
            do  
                if [[ "\$sample" ==  *"\$mod"* ]]
                then
                    echo "\$sample"
                    echo "\$sample" >> "\${output_location}/\${mod}_samples_vcf.txt"
                fi

                if [[ "\$mod" == "bulk_diss_polyA" && "\$sample" == *"2507"* ]]
                then
                    skip_sample=\$sample
                fi
            done
        done

        echo "Made lists"

        # Save variants separately for different modalities
        for mod in ${modalities.join(' ')}
        do
            echo "\$mod"
            bcftools view -S "\${output_location}/\${mod}_samples_vcf.txt" \
                -Oz -o "\${output_location}/\${mod}_modality_variants.vcf.gz" \
                $vcf_file
            echo "split"
            bcftools index \${output_location}/\${mod}_modality_variants.vcf.gz

            nr_samples=\$( bcftools query -l \${output_location}/\${mod}_modality_variants.vcf.gz | wc -l )
            read_depth_cutoff=\$(( ${params.read_depth} * nr_samples ))

            all_variants_output="\${output_location}/filtered_variants_rd_${params.read_depth}_\${mod}.vcf.gz"
            
            # Filter variants
            bcftools view \
                -Oz \
                -i "QUAL>=20 && DP>=\${read_depth_cutoff}" \
                -o \${all_variants_output} \
                \${output_location}/\${mod}_modality_variants.vcf.gz

            # Index the filtered VCF
            bcftools index \${all_variants_output}

            # Additional filter removing low quality sample 2507 from bulk_diss_polyA in hgsoc_new
            if [[ "\$mod" ==  "bulk_diss_polyA" ]]
            then
                bcftools view \
                    -Oz \
                    -s ^"\${skip_sample}" \
                    -o "\${output_location}/filtered_variants2_rd_${params.read_depth}_\${mod}.vcf.gz" \
                    \${all_variants_output}

                # Index the filtered VCF
                bcftools index "\${output_location}/filtered_variants2_rd_${params.read_depth}_\${mod}.vcf.gz"
            fi

        done
        
        echo "Finished splitting by modalities!"

        # Save variants individually
        for sample in `bcftools query -l $vcf_file` 
        do  
            echo "\$sample"
            sample_id="\${sample%????}"
            echo \$sample_id
            bcftools view \
                -c1 -Oz -s \$sample \
                -i "QUAL>=20 && DP>=${params.read_depth}" \
                -o "\${output_location}/\${sample_id}_individual_variants.vcf.gz" \
                $vcf_file

            bcftools index \${output_location}/\${sample_id}_individual_variants.vcf.gz
        done

        echo "Finished splitting output by sample"
        """
}