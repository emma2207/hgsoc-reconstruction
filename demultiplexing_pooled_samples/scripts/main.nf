#!/usr/bin/env nextflow

/*
 * Pipeline for demultiplexing pooled single-cell RNA-seq data
 */

nextflow.enable.dsl = 2

// Import processes
include { STAR_ALIGNMENT } from './modules/star_alignment.nf'
include { GENOTYPING } from './modules/genotyping.nf'
include { DEMULTIPLEXING } from './modules/demultiplexing.nf'


// Main workflow
workflow {
    // Input files
    fastq_pairs = channel.fromFilePairs(
        "${params.dataDir}/fastq/Pool*-GEX-*_R{1,2}.fastq.gz", 
        flat: true
        )
    barcode_files = channel.fromFilePairs(
        "${params.dataDir}/Pool*-GEX-*/outs/filtered_feature_bc_matrix/barcodes.tsv.gz", 
        flat: true
        )
    ref_vcf_file = channel.fromPath(
        "${params.refVcfDir}/${params.dataset}/real_data/ncells_null/read_depth_0/${params.ref_modality}_modality_variants.vcf.gz"
        )

    // Run STAR alignment
    star_input = fastq_pairs.map { sample_id, fastq_r1, fastq_r2 ->
            def matcher = (sample_id =~ /Pool(\d+)-.*/)
            assert matcher.matches(): "Could not parse pool from FASTQ sample id: ${sample_id}"

            def pool = matcher[0][1] as Integer
            tuple(pool, fastq_r1, fastq_r2)
        }
    aligned_reads = STAR_ALIGNMENT(star_input)

    // Run genotyping with cellSNP-lite
    genotyping_input = aligned_reads.map { pool, bam, bai ->
        tuple(pool, bam, bai, barcode_files.filter { i -> i.name.contains("Pool${pool}-") }.first())
    }
    genotyping_output = GENOTYPING(genotyping_input)

    // Run demultiplexing with Vireo
    demultiplexing_input = genotyping_output.map { pool, cellSNP_path ->
        tuple(pool, cellSNP_path, ref_vcf_file)
    }.set { demultiplexing_input }

    DEMULTIPLEXING(demultiplexing_input)

}