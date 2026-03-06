#!/usr/bin/env nextflow

process GENOTYPE_INDIVIDUAL {
    conda "${params.conda}/sample-matching"
    publishDir "${params.outdir}/1a_individual_vcf", mode: 'copy'
    
    input:
        path bam_file

    output:
        path("${params.dataset}/*/*/*.vcf.gz"), emit: vcf
        path("${params.dataset}/*/*/*.vcf.gz.csi"), emit: index

    script:
        """
        if [ ${params.pseudobulk} == true ]
        then 
            output_location="${params.dataset}/pseudobulk/ncells_${params.ncells}"
        else
            output_location="${params.dataset}/real_data/ncells_null"
        fi
        mkdir -p \$output_location
        
        sample_name=\$(basename ${bam_file} .bam)
        
        # Call variants using bcftools for individual sample
        bcftools mpileup -Ou -f ${params.refGenome}/fasta/genome.fa ${bam_file} | \
        bcftools call -mv -Ou | \
        bcftools view -Oz -i 'QUAL>=20' -o "\${output_location}/\${sample_name}.vcf.gz"

        # Index the VCF
        bcftools index "\${output_location}/\${sample_name}.vcf.gz"
        """
}
