#!/usr/bin/env nextflow

process PREPARE_INPUTS {
    tag "prepare inputs"
    conda "${params.conda}/sample-matching"
    publishDir "${params.outdir}/prepared_bams", mode: 'copy'
    
    input:
    tuple val(sample_id), val(nr), path(bam_files)
    
    output:
    path("*.bam"), emit: prepared_files

    script:
    """
    # Rename BAM files using the filename and the folder name
    # Original paths include "sample_id/pseudobulk_nr.bam"
    # The new names should be "sample_id_pseudobulk_nr.bam"
    for bam_file in ${bam_files}; do
        new_name="${sample_id}_${nr}.bam"
        mv "\$bam_file" "\$new_name"

        echo "Renamed \$bam_file to \$new_name"
    done
    
    # Temporary directory
    # tmpd=\$(mktemp -d)

    # Sort BAM files numerically and rename to "sample_id_pseudobulk_nr.bam"
    # for bam_file in ${bam_files}
    # do
    #    # extract header and split into HD / SQ / rest
    #    samtools view -H "\$bam_file" > "\$tmpd/header.sam"
    #    awk '
    #        /^@HD/   { print > "'"\$tmpd"'/hd.txt"; next }
    #        /^@SQ/   { print > "'"\$tmpd"'/sq.txt"; next }
    #        { print > "'"\$tmpd"'/rest.txt" }
    #    ' "\$tmpd/header.sam"

    #    # build new header: HD, sorted SQ (version/numeric aware), then other lines
    #    cat "\$tmpd/hd.txt" > "\$tmpd/new_header.sam"
    #    sort -t':' -k2,2 -V -f "\$tmpd/sq.txt" >> "\$tmpd/new_header.sam"
    #    cat "\$tmpd/rest.txt" >> "\$tmpd/new_header.sam"

    #    # write new BAM with reheader
    #    samtools reheader "\$tmpd/new_header.sam" "\$bam_file" > "${sample_id}_${nr}.bam"
    #done

    #rm -r "\$tmpd"
    """
}