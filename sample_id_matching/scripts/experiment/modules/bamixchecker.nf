nextflow.enable.dsl = 2

process BAMIXCHECKER {
    tag "${sample_id}"
    publishDir "${params.outdir}/bamixchecker", mode: 'copy'
    
    input:
      tuple val(sample_id), path(bam)
      path(ref)  // Reference genome fasta
      val(ref_ver) // Reference version (hg19 or hg38)
      path(bed) // Optional targeted regions BED file
    
    output:
      path("BAMixChecker/BAMixChecker_report.html"), emit: report
      path("BAMixChecker/BAMixChecker_heatmap.pdf"), emit: heatmap
      path("BAMixChecker/Total_result.txt"), emit: results
      path("BAMixChecker/Matched_samples.txt"), optional: true
      path("BAMixChecker/Mismatched_samples.txt"), optional: true

    script:
      def bed_arg = bed ? "--BEDfile ${bed}" : ""
      """
      set -euo pipefail

      # Create config file for BAMixChecker
      echo "GATK=gatk" > BAMixChecker.config
      echo "BEDTOOLS=bedtools" >> BAMixChecker.config

      # Create input file list
      echo "${bam}" > bam_list.txt

      # Run BAMixChecker
      python /path/to/BAMixChecker/BAMixChecker.py \\
        -l bam_list.txt \\
        -r ${ref} \\
        -v ${ref_ver} \\
        ${bed_arg} \\
        -o . \\
        -p 4
      """
}