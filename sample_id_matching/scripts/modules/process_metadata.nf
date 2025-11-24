#!/usr/bin/env nextflow

process PROCESS_METADATA {

      label 'process_metadata'

      publishDir 'results/'

      input:
      path(metadata)
      path(fastq_dir)

      output:
      path("run_name_mapping/run_name_mapping_${params.datatype}_${params.dataset}.csv")

      script:
      """
      python ${params.projectDir}/process_metadata.py \
            -d ${params.dataset} \
            -t ${params.datatype} \
            -m $metadata
      """
}