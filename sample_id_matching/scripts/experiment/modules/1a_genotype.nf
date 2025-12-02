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
        output_location="${params.dataset}"
        mkdir -p \$output_location

        # Call variants using bcftools and pipe directly to filtering
        bcftools mpileup -Ou -f ${params.refGenome}/fasta/genome.fa ${bam_files} | \
        bcftools call -mv -Ou | \
        bcftools view -Oz -i 'QUAL>=20 && DP>=1' -o \${output_location}/filtered_variants.vcf.gz

        # Index the filtered VCF
        bcftools index \${output_location}/filtered_variants.vcf.gz

        # Save variants individually
        for sample in `bcftools query -l \${output_location}/filtered_variants.vcf.gz` 
        do
            sample_id="\${sample%????}"
            echo \$sample_id
            bcftools view \
                -c1 -Oz -s \$sample \
                -o "\${output_location}/\${sample_id}_filtered_variants.vcf.gz" \
                \${output_location}/filtered_variants.vcf.gz

            bcftools index \${output_location}/\${sample_id}_filtered_variants.vcf.gz
        done
        """
}