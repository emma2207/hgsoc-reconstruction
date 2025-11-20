#!/usr/bin/env nextflow

process HYSYS {
    tag "hysys"
    conda "${params.conda}/hysys"
    publishDir "${params.outdir}/hysys", mode: 'copy'
    // errorStrategy 'ignore'
    
    input:
        path(vcf)
    
    output:
        path("concordance_output.txt"), emit: concordance
        path("model_results.txt"), emit: model
    
    script:
        """
        set -euo pipefail

        # Prep HYSYS input
        for file in ${vcf}
        do
            # get filename
            base_name=\$(basename "\$file")
            base_name=\${base_name%.vcf.gz}
            base_name=\${base_name%.vcf}

            # Add AF field to VCFs and convert VCFs to required format (SNP_ID CHR POS VAF)
            bcftools +fill-tags \$file -- -t AF | \
            bcftools query -f '%ID\t%CHROM\t%POS\t%AF\n' | \
            awk '\$4!="."' | sort -k1,1 > "\${base_name}.snps"

            # Add file-path to list
            echo "\${base_name}.snps" >> sample_list.txt
        done
        
        # HYSYS
        # Run concordance calculation (all vs all with symmetrical evaluation)
        ${params.HYSYS}/HaveYouSwappedYourSamples.sh conc -s sample_list.txt sample_list.txt concordance_output.txt
        
        # Run model analysis for detecting swaps
        ${params.HYSYS}/HaveYouSwappedYourSamples.sh model concordance_output.txt > model_results.txt
        """
}