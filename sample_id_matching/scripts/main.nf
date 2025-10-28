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
    metadata = channel.fromPath("${params.petaLibrary}/metadata/sra_run_tables/SraRunTable_${params.datatype}_${params.dataset}.csv")
    input_data = channel.fromPath("${params.petaLibrary}/${params.dataset}/${params.datatype}", type: 'dir')

    // Process metadata
    fastqs = PROCESS_METADATA(metadata, input_data)
    // fastqs.fastq_dirs_csv.view { x -> "Test channel: ${x}" }

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
            if (row.containsKey('R1_path_run_4')) {
                fastqs_r1 << file(row.'R1_path_run_4')
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
            if (row.containsKey('R2_path_run_4')) {
                fastqs_r2 << file(row.'R2_path_run_4')
            }

            tuple(row.'Sample Name', fastqs_r1, fastqs_r2)
        }
        .view { row -> "Sample fastq dirs: ${row}" }

    // QC
    fastp_out = QC_READS_WITH_FASTP(sample_ch)
    fastp_out.merged_fastqs.view { x -> "Fastp merged reads: ${x}" }
    
    // Align fastqs
    bam = ALIGNMENT_WITH_STAR(fastp_out.merged_fastqs)
    bam.aligned_reads.view { x -> "Aligned Reads: ${x}" }

    // Create pseudobulks
    n_barcodes = channel.of(params.n_barcodes)
    n_pseudobulks = channel.of(params.n_pseudobulks)
    pseudobulk_input = bam.aligned_reads.combine(n_barcodes).combine(n_pseudobulks)
    CREATE_PSEUDOBULKS(pseudobulk_input)
}
