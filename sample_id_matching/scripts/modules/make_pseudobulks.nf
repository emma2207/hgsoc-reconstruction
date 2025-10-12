#!/usr/bin/env nextflow

process CREATE_PSEUDOBULKS {
    
    label 'split'
    debug true

    input:
    tuple path(bam_dir), path(barcodes), val(n_barcodes), val(n_pseudobulks)

    output:
    true

    script:
    """
    mkdir -p "${params.outputDir}/subset-bam/${params.dataset}/${params.sample}"
    output_location="${params.outputDir}/subset-bam/${params.dataset}/${params.sample}"

    echo "Random barcode selection..."

    python ${params.projectDir}/select_barcodes.py \
        -d ${params.dataset} \
        -s ${params.sample} \
        -b "$bam_dir/Aligned.sortedByCoord.out.bam" \
        -n $n_barcodes \
        -j $n_pseudobulks

    for (( i=1; i<$n_pseudobulks+1; i++ ))
    do  
        echo "Subsetting pseudobulk \$i out of $n_pseudobulks"

        subset-bam -b "$bam_dir/Aligned.sortedByCoord.out.bam" \
            -c "${params.projectDir}/../output_data/subset-bam/${params.dataset}/${params.sample}/selected_barcodes_\$i.txt" \
            -o "pseudobulk_\$i.bam" \
            --cores 1 
    done

    echo "Finished subsetting!"

    cp *.bam \${output_location}/
    """

    stub:
    """
    echo "Splitting reads for sample ${params.sample} in dataset ${params.dataset}"
    """
}
