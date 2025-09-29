process QC_READS_WITH_FASTP {
      
      label 'fastp'

      input:
      path(fastq_r1)
      path(fastq_r2)

      output:
      tuple path("${params.sample}_R1_trimmed.fastq.gz"), 
            path("${params.sample}_R2_trimmed.fastq.gz")

      script:
      """
      mkdir -p "${params.outputDir}/fastp/${params.dataset}/${params.sample}"
      output_location="${params.outputDir}/fastp/${params.dataset}/${params.sample}"

      fastq_r1_array=($fastq_r1)
      fastq_r2_array=($fastq_r2)

      echo "Fastq files R1: \${fastq_r1_array[@]}"
      echo "Fastq files R2: \${fastq_r2_array[@]}"

      # Combine fastq files across lanes
      cat "\${fastq_r1_array[@]}" > "\${output_location}/${params.sample}_R1_merged.fastq.gz"
      cat "\${fastq_r2_array[@]}" > "\${output_location}/${params.sample}_R2_merged.fastq.gz"


      # Run fastp
      fastp --in1 "\${output_location}/${params.sample}_R1_merged.fastq.gz" \
            --in2 "\${output_location}/${params.sample}_R2_merged.fastq.gz" \
            --out1 "${params.sample}_R1_trimmed.fastq.gz" \
            --out2 "${params.sample}_R2_trimmed.fastq.gz" \
            --html "\${output_location}/${params.sample}_fastp_report.html" \
            --json "\${output_location}/${params.sample}_fastp_report.json" \
            --thread 3
      """

      stub:
      """
      echo "Running fastp on ${fastq_r1} and ${fastq_r2} for sample ${params.sample} \
      in dataset ${params.dataset}"
      """
}
