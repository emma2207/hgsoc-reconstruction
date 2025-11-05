nextflow.enable.dsl = 2

process FILTER_BAM {
    tag "${sample_id}"

    publishDir "${params.outdir}/filtered_bams", mode: 'copy'

    input:
      tuple val(sample_id), path(bam)
    output:
      tuple val(sample_id), path("${sample_id}.filtered.bam"), path("${sample_id}.filtered.bam.bai")
    
    script:
      """
      set -euo pipefail
      
      # Filter BAM file for:
      # -q 20: minimum mapping quality of 20
      # -F 3844: remove unmapped, not primary alignment, supplementary alignments, duplicates
      # -@ 4: use 4 threads
      samtools view -bh -q 20 -F 3844 -@ 4 ${bam} > ${sample_id}.filtered.bam

      # Index the filtered BAM
      samtools index ${sample_id}.filtered.bam
      """
}