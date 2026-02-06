#!/usr/bin/env nextflow

process GENOTYPE_AND_FILTER {
    conda "${params.conda}/sample-matching"
    publishDir "${params.outdir}/1a_vcf", mode: 'copy'
    
    input:
        path bam_files
        val modalities 

    
    output:
        path("${params.dataset}/*/*/*/*_modality_variants.vcf.gz"), emit: modality_vcfs
        path("${params.dataset}/*/*/*/*_individual_variants.vcf.gz"), emit: individual_vcfs
        path("${params.dataset}/*/*/*/filtered_variants.vcf.gz")
        path("${params.dataset}/*/*/*/filtered_variants_rd_*.vcf.gz")

    script:
        """
        if [ ${params.pseudobulk} == true ]
        then 
            output_location="${params.dataset}/pseudobulk/ncells_${params.ncells}/read_depth_1"
        else
            output_location="${params.dataset}/real_data/ncells_null/read_depth_1" 
        fi
        all_variants_output="\${output_location}/filtered_variants.vcf.gz"
        rd_filtered_variants="\${output_location}/filtered_variants_rd_1.vcf.gz"
        mkdir -p \$output_location

        # List samples by modality
        for mod in ${modalities.join(' ')}
        do  
            echo "\$mod"
            for bam in $bam_files
            do
                echo "\$bam"
                if [[ "\$bam" == *"\$mod"* ]]
                then
                    echo "\$bam" >> "\${output_location}/\${mod}_files.txt"
                fi
            done
        done

        echo "Start genotyping"

        # Call variants using bcftools and pipe directly to filtering
        bcftools mpileup -Ou -f ${params.refGenome}/fasta/genome.fa ${bam_files} | \
        bcftools call -mv -Ou | \
        bcftools view -Oz -i 'QUAL>=20' -o \${all_variants_output}

        echo "Finished genotyping!"

        # Index the filtered VCF
        bcftools index \${all_variants_output}

        echo "Finished indexing..."

        # Filter variants further
        bcftools view -Oz -i 'DP>=1' -o \${rd_filtered_variants} \${all_variants_output}
        bcftools index \${rd_filtered_variants}

        echo "Filtered variants by read depth"

        # Save variants separately for different modalities
        for mod in ${modalities.join(' ')}
        do
            bcftools view -S "\${output_location}/\${mod}_files.txt" \
                -Oz -o "\${output_location}/\${mod}_modality_variants.vcf.gz" \
                \${rd_filtered_variants}
            bcftools index \${output_location}/\${mod}_modality_variants.vcf.gz
        done

        echo "Saved variants per modality"

        # Save variants individually
        for sample in `bcftools query -l \${rd_filtered_variants}`  
        do
            sample_id="\${sample%????}"
            echo \$sample_id
            bcftools view \
                -c1 -Oz -s \$sample \
                -o "\${output_location}/\${sample_id}_individual_variants.vcf.gz" \
                \${rd_filtered_variants}

            bcftools index \${output_location}/\${sample_id}_individual_variants.vcf.gz
        done

        echo "Finished splitting output by sample"
        """
}