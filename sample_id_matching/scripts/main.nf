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
    // Input metadata & fastqs
    metadata = Channel.fromPath("${params.dataDir}/metadata/SraRunTable_${params.datatype}_${params.dataset}_test.csv")
    input_data = Channel.fromPath("${params.dataDir}/test_fastq/${params.dataset}/${params.datatype}", type: 'dir')

    // Process metadata
    fastqs = PROCESS_METADATA(metadata, input_data)
    fastqs.fastq_dirs_csv.view { x -> "Test channel: ${x}" }

    // Split the fastqs by R1 and R2 paths for each sample
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

            tuple(row.'Sample Name', fastqs_r1, fastqs_r2)
        }
        .view { row -> "Sample fastq dirs: ${row}" }

    sample_name = sample_ch.map { it[0] }
    sample_name.view { x -> "Sample name: ${x}" }

    // I need to extract the sample names from the csv and feed those to the subsequent processes

    // QC
    fastp_out = QC_READS_WITH_FASTP(sample_ch)
    fastp_out.merged_R1.view { x -> "Fastp merged R1: ${x}" }
    
    // Align fastqs
    merged_fastqs = fastp_out.merged_R1.combine(fastp_out.merged_R2)
    merged_fastqs.view { x -> "Merged fastqs: ${x}" }
    bam = ALIGNMENT_WITH_STAR(sample_name.combine(merged_fastqs))
    bam.aligned_reads.view { x -> "Aligned Reads: ${x}" }

    // Create pseudobulks
    n_barcodes = channel.of(100)
    n_pseudobulks = channel.of(1)
    CREATE_PSEUDOBULKS(sample_name.combine(bam.aligned_reads), n_barcodes, n_pseudobulks)
}
