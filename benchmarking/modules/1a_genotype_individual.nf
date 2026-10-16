#!/usr/bin/env nextflow

process GENOTYPE_INDIVIDUAL {
    conda "${params.conda}/sample-matching"
    publishDir "${params.outdir}/1a_individual_vcf", mode: 'copy'
    tag "${bam_file.simpleName}"
    
    input:
        tuple path(bam_file), path(bai_file)
        path(snp_ref_vcf)
        path(snp_ref_tbi)

    output:
        path("${params.dataset}/*/*/*.vcf.gz"), emit: vcf
        path("${params.dataset}/*/*/*.vcf.gz.csi"), emit: index

    script:
        """
        set -euo pipefail

        echo "START" >&2
        date -Is >&2
        hostname >&2
        env | grep '^SLURM_JOB_ID=' >&2 || true
        echo "BAM=${bam_file}" >&2
        df -h . >&2 || true

        if [ "${params.pseudobulk}" = "true" ]
        then 
            output_location="${params.dataset}/pseudobulk/ncells_${params.ncells}"
        else
            output_location="${params.dataset}/real_data/ncells_null"
        fi
        mkdir -p \$output_location
        
        sample_name=\$(basename ${bam_file} .bam)
        
        # Call variants using bcftools for individual sample
        bcftools mpileup -Ou \
            -d ${params.mpileup_max_depth} \
            -f ${params.refGenome}/fasta/genome.fa \
            -T ${snp_ref_vcf} \
            ${bam_file} | \
        bcftools call -m -Ou | \
        bcftools view -Oz -i 'QUAL>=20' -o "\${output_location}/\${sample_name}.vcf.gz"

        # Index only when a non-empty VCF was produced.
        if [ -s "\${output_location}/\${sample_name}.vcf.gz" ]; then
            bcftools index -f "\${output_location}/\${sample_name}.vcf.gz"
        else
            echo "Skipping index: \${output_location}/\${sample_name}.vcf.gz is empty or missing"
        fi

        echo "DONE" >&2
        date -Is >&2
        """
}
