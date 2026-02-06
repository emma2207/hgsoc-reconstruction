#!/usr/bin/env nextflow

process PREPARE_INPUTS {
    conda "${params.conda}/sample-matching"
    publishDir "${params.outdir}/0_prepared_bams", mode: 'copy'
    
    input:
        tuple val(sample_id), val(nr), path(bam_files), val(datatype)
    
    output:
        path("${params.dataset}/*/*/*.bam")

    script:
        """
        if [ ${params.pseudobulk} == true ]
        then 
            output_location="${params.dataset}/pseudobulk/ncells_${params.ncells}"
        else
            output_location="${params.dataset}/real_data/ncells_null" 
        fi
        mkdir -p \$output_location
        
        # Rename BAM files using the filename and the folder name
        if [ ${params.pseudobulk} == true ]
        then
            # Original pseudobulk paths are "sample_id/pseudobulk_nr.bam"
            # New names are "sample_id_pseudobulk_nr.bam"
            for bam_file in ${bam_files}
            do
                new_name="\${output_location}/${sample_id}_${nr}.bam"
                mv "\$bam_file" "\$new_name"

                echo "Renamed \$bam_file to \$new_name"
            done
        else
            # Original bulk paths are "datatype/sample_id/Aligned.sortedByCoord.out.bam"
            # New names are "sample_id_datatype.bam"
            for bam_file in ${bam_files}
            do
                new_name="\${output_location}/${sample_id}_${datatype}.bam"
                mv "\$bam_file" "\$new_name"

                echo "Renamed \$bam_file to \$new_name"
            done
        fi
        """
}