process QC_READS_WITH_FASTP {
      
      label 'fastp'

      publishDir 'results/', mode: 'copy'
      
      errorStrategy 'ignore'

      input:
      tuple val(sample_name), path(R1_paths), path(R2_paths)

      output:
      path("fastp/${params.dataset}/${params.datatype}/${sample_name}_fastp_report.html")
      path("fastp/${params.dataset}/${params.datatype}/${sample_name}_fastp_report.json")
      tuple val(sample_name),
            path("fastp/${params.dataset}/${params.datatype}/${sample_name}_R1_merged.fastq.gz"),
            path("fastp/${params.dataset}/${params.datatype}/${sample_name}_R2_merged.fastq.gz"), emit: merged_fastqs

      script:
      """
      output_location="fastp/${params.dataset}/${params.datatype}"
      mkdir -p \$output_location

      echo "Fastq files R1: $R1_paths"
      echo "Fastq files R2: $R2_paths"

      # Combine fastq files across lanes
      cat $R1_paths > "\${output_location}/${sample_name}_R1_merged.fastq.gz"
      cat $R2_paths > "\${output_location}/${sample_name}_R2_merged.fastq.gz"

      # Run fastp
      fastp --in1 "\${output_location}/${sample_name}_R1_merged.fastq.gz" \
            --in2 "\${output_location}/${sample_name}_R2_merged.fastq.gz" \
            --html "\${output_location}/${sample_name}_fastp_report.html" \
            --json "\${output_location}/${sample_name}_fastp_report.json" \
            --disable_adapter_trimming \
            --thread 3
      """
}
