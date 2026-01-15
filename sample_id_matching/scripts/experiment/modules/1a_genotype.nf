#!/usr/bin/env nextflow

process GENOTYPE_AND_FILTER {
    conda "${params.conda}/sample-matching"
    publishDir "${params.outdir}/vcf", mode: 'copy'
    
    input:
        path bam_files
        val modalities 

    
    output:
        // path("${params.dataset}/*/filtered_variants.vcf.gz"), emit: combined_vcf
        path("${params.dataset}/*/filtered_variants_*.vcf.gz"), emit: modality_vcfs
        path("${params.dataset}/*/*_filtered_variants.vcf.gz"), emit: individual_vcfs

    script:
        """
        if [ ${params.pseudobulk} == true ]
        then 
            output_location="${params.dataset}/pseudobulk"
            all_variants_output="\${output_location}/filtered_variants.vcf.gz"
            rd_filtered_variants="\${output_location}/rd_filtered_variants.vcf.gz"
        else
            output_location="${params.dataset}/real_data"
            all_variants_output="\${output_location}/filtered_variants.vcf.gz"
            rd_filtered_variants="\${output_location}/rd_filtered_variants.vcf.gz"
        fi
        output_location="${params.dataset}/real_data"
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

        # Save variants separately for different modalities
        for mod in ${modalities.join(' ')}
        do
            bcftools view -S "\${output_location}/\${mod}_files.txt" \
                -Oz -o "\${output_location}/filtered_variants_\${mod}.vcf.gz" \
                \${rd_filtered_variants}
            bcftools index \${output_location}/filtered_variants_\${mod}.vcf.gz
        done

        # Save variants individually
        for sample in `bcftools query -l \${rd_filtered_variants}`  
        do
            sample_id="\${sample%????}"
            echo \$sample_id
            bcftools view \
                -c1 -Oz -s \$sample \
                -o "\${output_location}/\${sample_id}_filtered_variants.vcf.gz" \
                \${rd_filtered_variants}

            bcftools index \${output_location}/\${sample_id}_filtered_variants.vcf.gz
        done

        echo "Finished splitting output by sample"
        """
}