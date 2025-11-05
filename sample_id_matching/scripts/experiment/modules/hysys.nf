nextflow.enable.dsl = 2

process HYSYS {
    tag "${sample_id}"
    publishDir "${params.outdir}/hysys", mode: 'copy'
    
    input:
      tuple val(sample_id), path(vcf)
    
    output:
      path("concordance_output.txt"), emit: concordance
      path("hysys_heatmap.pdf"), emit: heatmap
      path("model_results.txt"), emit: model
    
    script:
      """
      set -euo pipefail
      
      # Convert VCF to required format (SNP_ID CHR POS VAF)
      bcftools query -f '%ID\t%CHROM\t%POS\t%AF\n' ${vcf} | \
      awk '\$4!="."' | sort -k1,1 > ${sample_id}.snps
      
      # Create sample list file
      echo "\$PWD/${sample_id}.snps" > sample_list.txt
      
      # Run concordance calculation (all vs all with symmetrical evaluation)
      HaveYouSwappedYourSamples.sh conc -s sample_list.txt sample_list.txt concordance_output.txt
      
      # Generate heatmap
      HaveYouSwappedYourSamples.sh heat concordance_output.txt hysys_heatmap.pdf
      
      # Run model analysis for detecting swaps
      HaveYouSwappedYourSamples.sh model concordance_output.txt > model_results.txt
      """
}