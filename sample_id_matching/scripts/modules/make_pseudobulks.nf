#!/usr/bin/env nextflow

process CREATE_PSEUDOBULKS {
    
    label 'pseudobulks'

    publishDir 'results/', mode: 'copy'

    input:
    tuple val(sample_name), 
        path(bam), 
        path(bai),
        val(n_barcodes),
        val(n_pseudobulks)

    output:
    path("subset-bam/${params.dataset}/${sample_name}/pseudobulk_*.bam")

    script:
    """
    output_location="subset-bam/${params.dataset}/${sample_name}"
    mkdir -p \$output_location

    echo "Random barcode selection..."

    python ${params.projectDir}/select_barcodes.py \
        -d ${params.dataset} \
        -s ${sample_name} \
        -b "$bam" \
        -n $n_barcodes \
        -j $n_pseudobulks

    for (( i=1; i<$n_pseudobulks+1; i++ ))
    do  
        echo "Subsetting pseudobulk \$i out of $n_pseudobulks"

        subset-bam -b "$bam" \
            -c "subset-bam/${params.dataset}/${sample_name}/selected_barcodes_\$i.txt" \
            -o "\${output_location}/pseudobulk_\$i.bam" \
            --cores 1 
    done

    echo "Finished subsetting!"

    """
}
