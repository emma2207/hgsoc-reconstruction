#!/usr/bin/env nextflow

process GENOTYPE_AND_FILTER {
    conda "${params.conda}/sample-matching"
    publishDir "${params.outdir}/vcf", mode: 'copy'
    
    input:
        path(bam_files)

    
    output:
        path("${params.dataset}/filtered_variants.vcf.gz"), emit: combined_vcf
        path("${params.dataset}/*_filtered_variants.vcf.gz"), emit: individual_vcfs

    script:
        """
        if [ ${params.pseudobulk} == true ]
        then 
            output_location="${params.dataset}/pseudobulk"
            all_variants_output="\${output_location}/filtered_variants.vcf.gz"
        else
            output_location="${params.dataset}/real_data"
            all_variants_output="\${output_location}/filtered_variants.vcf.gz"
        fi
        mkdir -p \$output_location

        echo "Start genotyping"

        # Call variants using bcftools and pipe directly to filtering
        bcftools mpileup -Ou -f ${params.refGenome}/fasta/genome.fa ${bam_files} | \
        bcftools call -mv -Ou | \
        bcftools view -Oz -i 'QUAL>=20 && DP>=1' -o \${all_variants_output}

        echo "Finished genotyping!"

        # Index the filtered VCF
        bcftools index \${all_variants_output}

        echo "Finished indexing..."

        # Save variants individually
        for sample in `bcftools query -l \${all_variants_output}` 
        do
            sample_id="\${sample%????}"
            echo \$sample_id
            bcftools view \
                -c1 -Oz -s \$sample \
                -o "\${output_location}/\${sample_id}_filtered_variants.vcf.gz" \
                \${all_variants_output}

            bcftools index \${output_location}/\${sample_id}_filtered_variants.vcf.gz
        done

        echo "Finished splitting output by sample"
        """
}