#!/usr/bin/env nextflow

/*
 * Pipeline for demultiplexing pooled single-cell RNA-seq data
 */

nextflow.enable.dsl = 2

def extractPool(value, label) {
    def matcher = (value.toString() =~ /(?:^|\/)Pool(\d+)-/)
    assert matcher.find(): "Could not parse pool from ${label}: ${value}"
    matcher.group(1) as Integer
}

// Import processes
include { STAR_ALIGNMENT } from './modules/star_alignment.nf'
include { GENOTYPING } from './modules/genotyping.nf'
include { DEMULTIPLEXING } from './modules/demultiplexing.nf'
include { SPLIT_POOLED_BAM } from './modules/split_pooled_bam.nf'


// Main workflow
workflow {
    // Input files
    fastq_pairs = channel.fromFilePairs(
        "${params.dataDir}/fastq/Pool*-GEX-*_R{1,2}_*.fastq.gz", 
        flat: true
        )
    barcode_files = channel.fromPath(
        "${params.dataDir}/Pool*-GEX-*/outs/filtered_feature_bc_matrix/barcodes.tsv.gz", 
        ).map { barcode_file ->
            def pool = extractPool(barcode_file, 'barcode path')
            tuple(pool, barcode_file)
        }
    ref_vcf_file = channel.fromPath(
        "${params.refVcfDir}/${params.dataset}/real_data/ncells_null/read_depth_0/${params.ref_modality}_modality_variants.vcf.gz"
        )

    // Run STAR alignment
    star_input = fastq_pairs.map { sample_id, fastq_r1, fastq_r2 ->
            def pool = extractPool(sample_id, 'FASTQ sample id')
            tuple(pool, fastq_r1, fastq_r2)
        }
        .groupTuple()
        .map { pool, fastq_r1_list, fastq_r2_list ->
            def sorted_r1 = fastq_r1_list.sort { fastq -> fastq.name }
            def sorted_r2 = fastq_r2_list.sort { fastq -> fastq.name }

            assert sorted_r1.size() == sorted_r2.size(): "Mismatched R1/R2 lane counts for pool ${pool}"
            tuple(pool, sorted_r1, sorted_r2)
        }
    star_input.view { pool, r1, r2 -> println "Pool ${pool}: ${r1.size()} R1 files, ${r2.size()} R2 files" }
    aligned_reads = STAR_ALIGNMENT(star_input)

    // Run genotyping with cellSNP-lite
    genotyping_input = aligned_reads.join(barcode_files)
        .map { pool, bam, bai, barcodes -> tuple(pool, bam, bai, barcodes) }
    genotyping_output = GENOTYPING(genotyping_input)

    // Run demultiplexing with Vireo
    demultiplexing_input = genotyping_output
        .map { pool, cellSNP_path -> tuple(pool, cellSNP_path) }
    demultiplexing_output = DEMULTIPLEXING(demultiplexing_input, ref_vcf_file)

    // Split pooled BAM into per-donor BAM files after demultiplexing
    split_bam_input = aligned_reads.join(demultiplexing_output)
        .map { pool, bam, bai, vireo_pool_dir -> tuple(pool, bam, bai, vireo_pool_dir) }
    SPLIT_POOLED_BAM(split_bam_input)

}