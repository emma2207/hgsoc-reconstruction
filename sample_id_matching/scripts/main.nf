#!/usr/bin/env nextflow

/*
 * main.nf — run simulated data pipeline
 * 
 * Every process has to be run for every sample in every dataset that we want to create
 * pseudobulks for.
 */

nextflow.enable.dsl = 2

include { QC_READS_WITH_FASTP } from './processes/fastp.nf'
// include { ALIGNMENT_WITH_STAR } from './processes/star_alignment.nf'
// include { SPLIT_READS } from './processes/subset-bam.nf'
// include { CREATE_PSEUDOBULK } from './processes/pseudobulk_creation.nf'

workflow {
    // Create initial trigger channel
    // fastq_r1 = Channel.fromPath("${params.petaLibrary}/Pool1-GEX-*_R1_001.fastq.gz").collect()
    // fastq_r2 = Channel.fromPath("${params.petaLibrary}/Pool1-GEX-*_R2_001.fastq.gz").collect()
    fastq_r1 = Channel.fromPath("${params.testSet}/${params.sample}/SRR*_1_truncated.fastq.gz").collect()
    fastq_r2 = Channel.fromPath("${params.testSet}/${params.sample}/SRR*_2_truncated.fastq.gz").collect()

    // fastqs = fastq_r1.combine(fastq_r2).collect()
    // Run pipeline processes
    QC_READS_WITH_FASTP(fastq_r1, fastq_r2)
    // ALIGNMENT_WITH_STAR(ref_genome, QC_READS_WITH_FASTP.out)
    // SPLIT_READS(cell_barcodes, ALIGNMENT_WITH_STAR.out)
    // CREATE_PSEUDOBULK(SPLIT_READS.out)
}
