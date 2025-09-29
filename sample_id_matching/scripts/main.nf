#!/usr/bin/env nextflow

/*
 * main.nf — run simulated data pipeline
 */

nextflow.enable.dsl = 2

include { QC_READS_WITH_FASTP } from './modules/fastp.nf'
include { ALIGNMENT_WITH_STAR } from './modules/star_alignment.nf'
include { CREATE_PSEUDOBULKS } from './modules/make_pseudobulks.nf'

workflow {
    // Input fastqs
    fastq_r1 = Channel.fromPath("${params.testSet}/${params.sample}/SRR*_1_truncated.fastq.gz").collect()
    fastq_r2 = Channel.fromPath("${params.testSet}/${params.sample}/SRR*_2_truncated.fastq.gz").collect()

    // Input barcodes
    barcodes = Channel.fromPath("${params.projectDir}/../data/barcodes/${params.dataset}/${params.sample}/*_sub.csv", checkIfExists: true)

    // QC
    trimmed_fastqs = QC_READS_WITH_FASTP(fastq_r1, fastq_r2)
    // Input merged fastqs
    fastq_r1 = Channel.fromPath("${params.outputDir}/fastp/${params.dataset}/${params.sample}/${params.sample}_R1_merged.fastq.gz").collect()
    fastq_r2 = Channel.fromPath("${params.outputDir}/fastp/${params.dataset}/${params.sample}/${params.sample}_R2_merged.fastq.gz").collect()
    merged_fastqs = fastq_r1.combine(fastq_r2)
    // Align fastqs
    bam = ALIGNMENT_WITH_STAR(merged_fastqs)

    bam_dir = Channel.fromPath("${params.outputDir}/star/${params.dataset}/${params.sample}", type: 'dir')
    n_barcodes = channel.of(100)
    n_pseudobulks = channel.of(1)
    split_reads_input = bam_dir.combine(barcodes).combine(n_barcodes).combine(n_pseudobulks) 
    split_reads_input.view()
    // Create pseudobulks
    barcodes = CREATE_PSEUDOBULKS(split_reads_input)
}
