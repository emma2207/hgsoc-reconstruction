#!/usr/bin/env nextflow

process GENOTYPE_AND_FILTER {
    tag "${sample_id}"
    publishDir "${params.outdir}/vcfs", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam_file)
    
    output:
    tuple val(sample_id), path("${sample_id}.filtered.vcf"), emit: filtered_vcf
    
    script:
    """
    # Call variants using bcftools and pipe directly to filtering
    bcftools mpileup -Ou -f ${params.genome} ${bam_file} | \
    bcftools call -mv -Ou | \
    bcftools filter -i 'QUAL>=20 && DP>=10' -Ov -o ${sample_id}.filtered.vcf

    # Index the filtered VCF
    bcftools index ${sample_id}.filtered.vcf
    """

    stub:
    """
    touch ${sample_id}.filtered.vcf
    """
}