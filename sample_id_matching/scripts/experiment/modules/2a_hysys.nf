#!/usr/bin/env nextflow

process HYSYS {
    tag "hysys"
    conda "${params.conda}/hysys"
    publishDir "${params.outdir}/hysys", mode: 'copy'
    errorStrategy 'ignore'
    
    input:
        path(vcf)
    
    output:
        path("${params.dataset}/concordance_output.txt")
        path("${params.dataset}/model_results.txt")
    
    script:
        """
        set -euo pipefail

        output_location="${params.dataset}"
        mkdir -p \$output_location

        # Prep HYSYS input
        for file in ${vcf}
        do
            # get filename
            base_name=\$(basename "\$file")
            base_name=\${base_name%.vcf.gz}
            base_name=\${base_name%.vcf}

            # Add AF field to VCFs and convert VCFs to required format 
            # (SNP_ID CHR POS VAF)
            bcftools +fill-tags \$file -- -t AF | \
            bcftools query -f '%ID\t%CHROM\t%POS\t%AF\n' | \
            awk '\$4!="."' | sort -k1,1 > "\${base_name}.snps"

            # Add file-path to list
            echo "\${base_name}.snps" >> \${output_location}/sample_list.txt
        done
        
        # HYSYS
        # Run concordance calculation
        ${params.HYSYS}/HaveYouSwappedYourSamples.sh conc -s \\
            \${output_location}/sample_list.txt \\
            \${output_location}/sample_list.txt \\
            \${output_location}/concordance_output.txt
            
        # Run model analysis for detecting swaps
        ${params.HYSYS}/HaveYouSwappedYourSamples.sh model \\
            \${output_location}/concordance_output.txt \\
            > \${output_location}/model_results.txt
        """
}