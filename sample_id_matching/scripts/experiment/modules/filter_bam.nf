#!/usr/bin/env nextflow

process FILTER_BAM {
    tag "filter bam"
    conda "${params.conda}/sample-matching"
    publishDir "${params.outdir}/filtered_bam", mode: 'copy'

    input:
      tuple val(sample_id), val(pseudobulk), path(bam_files)
    
    output:
        path("${params.dataset}/${sample_id}_${pseudobulk}_filtered.bam"), emit: bam
        path("${params.dataset}/${sample_id}_${pseudobulk}_filtered.bam.bai")
    
    script:
      """
      set -euo pipefail

      output_location="${params.dataset}"
      mkdir -p \$output_location
      
      # Filter BAM file
      # -F 3844: remove unmapped, not primary alignment, supplementary alignments, duplicates
      samtools view -bh -q 20 -F 3844 ${bam_files} > \${output_location}/${sample_id}_${pseudobulk}_filtered.bam

      # Index the filtered BAM
      samtools index \${output_location}/${sample_id}_${pseudobulk}_filtered.bam
      """
}