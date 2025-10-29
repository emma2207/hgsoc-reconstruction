#!/usr/bin/env nextflow

process GENOTYPE_AND_FILTER {
    tag "${sample_id}"
    publishDir "${params.outdir}/bcftools", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam_file)
    
    output:
    tuple val(sample_id), path("${sample_id}_filtered.vcf"), emit: filtered_vcf

    script:
    """
    # Call variants using bcftools and pipe directly to filtering
    bcftools mpileup -Ou -f ${params.genome}/fasta/genome.fa ${bam_file} | \
    bcftools call -mv -Ov | \
    bcftools view -i 'QUAL>=20 && DP>=30' -o ${sample_id}_filtered.vcf

    # Index the filtered VCF
    bcftools index ${sample_id}_filtered.vcf
    """
}