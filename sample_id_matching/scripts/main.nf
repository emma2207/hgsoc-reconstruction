#!/usr/bin/env nextflow

/*
 * main.nf — run simulated data pipeline
 */

nextflow.enable.dsl = 2

include { PROCESS_METADATA } from './modules/process_metadata.nf'
include { QC_READS_WITH_FASTP } from './modules/fastp.nf'
include { ALIGNMENT_WITH_STAR } from './modules/star_alignment.nf'
include { CREATE_PSEUDOBULKS } from './modules/make_pseudobulks.nf'

workflow {
    // Input metadata
    metadata = Channel.fromPath("${params.dataDir}/metadata/SraRunTable_${params.datatype}_${params.dataset}_test.csv")
    input_data = Channel.fromPath("${params.dataDir}/test_fastq/${params.dataset}/${params.datatype}", type: 'dir')

    // Process metadata
    fastqs = PROCESS_METADATA(metadata, input_data)
    // fastqs.fastq_dirs_csv.view { x -> "Test channel: ${x}" }

    // Note that the fastqs directory contains an arbitrary number of samples, all of which we want to process 
    // through the rest of the pipeline on their own.
    sample_ch = fastqs.fastq_dirs_csv
        .splitCsv( header: true )
        .map { row ->
            def fastqs_r1 = []
            def fastqs_r2 = []
            if (row.containsKey('R1_path_run_0')) {
                fastqs_r1 << file(row.'R1_path_run_0')
            } 
            if (row.containsKey('R1_path_run_1')) {
                fastqs_r1 << file(row.'R1_path_run_1')
            } 
            if (row.containsKey('R1_path_run_2')) {
                fastqs_r1 << file(row.'R1_path_run_2')
            } 
            if (row.containsKey('R1_path_run_3')) {
                fastqs_r1 << file(row.'R1_path_run_3')
            } 
            if (row.containsKey('R2_path_run_0')) {
                fastqs_r2 << file(row.'R2_path_run_0')
            }
            if (row.containsKey('R2_path_run_1')) {
                fastqs_r2 << file(row.'R2_path_run_1')
            }
            if (row.containsKey('R2_path_run_2')) {
                fastqs_r2 << file(row.'R2_path_run_2')
            }
            if (row.containsKey('R2_path_run_3')) {
                fastqs_r2 << file(row.'R2_path_run_3')
            }

            tuple(fastqs_r1, fastqs_r2)
        }
        .view { row -> "Sample fastq dirs: ${row}" }

    // QC
    fastp_out = QC_READS_WITH_FASTP(sample_ch)
    fastp_out.merged_R1.view { x -> "Fastp merged R1: ${x}" }

    // Input barcodes
    // barcodes = Channel.fromPath("${params.projectDir}/../data/barcodes/${params.dataset}/${params.sample}/*_sub.csv", checkIfExists: true)

    
    // Input merged fastqs
    merged_fastqs = fastp_out.merged_R1.combine(fastp_out.merged_R2)
    merged_fastqs.view { x -> "Merged fastqs: ${x}" }
    // Align fastqs
    bam = ALIGNMENT_WITH_STAR(merged_fastqs)

    // bam_dir = Channel.fromPath("${params.outputDir}/star/${params.dataset}/${params.sample}", type: 'dir')
    
    bam.aligned_reads.view { x -> "Aligned Reads: ${x}" }
    n_barcodes = channel.of(100)
    n_pseudobulks = channel.of(1)
    // // Create pseudobulks
    barcodes = CREATE_PSEUDOBULKS(bam.aligned_reads, n_barcodes, n_pseudobulks)
}
