#!/usr/bin/env nextflow

process FILTER_BAM {
    conda "${params.conda}/sample-matching"
    publishDir "${params.outdir}/1b_filtered_bam", mode: 'copy'

    input:
      tuple val(sample_id), val(pseudobulk), path(bam_files), val(datatype)
    
    output:
        path("${params.dataset}/*/*_filtered.bam"), emit: bam

    script:
        """
        set -euo pipefail

        if [ ${params.pseudobulk} == true ]
        then 
            output_location="${params.dataset}/pseudobulk/ncells_${params.ncells}"
            output_name=\${output_location}/${sample_id}_${pseudobulk}_filtered.bam
        else
            output_location="${params.dataset}/real_data"
            output_name=\${output_location}/${sample_id}_${datatype}_filtered.bam
        fi
        mkdir -p \$output_location
        
        # Filter BAM file
        # -F 3844: remove unmapped, not primary alignment, supplementary alignments, duplicates
        samtools view -bh -q 20 -F 3844 ${bam_files} > \${output_name}

        # Index the filtered BAM
        samtools index \${output_name}
        """
}