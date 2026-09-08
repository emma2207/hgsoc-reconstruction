#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PROCESS_METADATA } from './modules/process_metadata.nf'
include { PSEUDOBULK_BAM_TO_FASTQ } from './modules/pseudobulk_bam_to_fastq.nf'
include { NTSM_COUNT } from './modules/ntsm_count.nf'
include { NTSM_EVAL } from './modules/ntsm_eval.nf'

// Collect all columns matching a run prefix (e.g., R1_path_run_0, R1_path_run_1)
// and return them in numeric run order.
def collectRunsDynamic(row, prefix) {
    // collect keys and build a plain-string regex (no slashy regex interpolation)
    def keys   = row.keySet().toList()
    def regex  = "${prefix}_run_\\d+"        // note: double-escape for \d inside a GString
    def matched = keys.findAll { k -> k ==~ regex }
    def sorted = matched.sort { a, b ->
        (a =~ /\d+$/)[0].toInteger() <=> (b =~ /\d+$/)[0].toInteger()
    }
    def values = sorted.collect { k -> row[k]?.toString()?.trim() }.findAll { value -> value }
    def files = values.collect { value -> file(value) }

    return files
}

workflow {
    if (!params.pseudobulk) {
        modalities = ["bulk_dissociated_polyA", "bulk_chunk_ribo"]
        metadata_patterns = modalities.collect { mod ->
            "${params.petaLibrary}/metadata/sra_run_tables/SraRunTable_${mod}_${params.dataset}_${params.read_type}.csv"
        }
        fastq_patterns = modalities.collect { mod ->
            "${params.petaLibrary}/raw_reads/${params.dataset}/${mod}"
        }
        metadata_ch = channel.fromPath(metadata_patterns)
        fastq_dir_ch = channel.fromPath(fastq_patterns, type: 'dir')
        metadata_ch.view { x -> "Metadata files: ${x}"}
        fastq_dir_ch.view { x -> "Fastq files: ${x}" }

        mapping_csv = PROCESS_METADATA(metadata_ch, fastq_dir_ch)

        sample_fastqs = mapping_csv.splitCsv(header: true)
                .map { row ->
                    tuple(
                        row.'Sample Name',
                        collectRunsDynamic(row, 'R1_path'),
                        collectRunsDynamic(row, 'R2_path')
                    )
                }
            .view { row -> "Sample fastq dirs: ${row}" }
    } 
    else {
        modalities = ["pseudobulk"]
        pseudobulk_bam_pattern = "${params.prepared_bam_dir}/${params.dataset}/pseudobulk/ncells_${params.ncells}/*.bam"
        pseudobulk_bams = channel
            .fromPath(pseudobulk_bam_pattern)
            .map { bam -> tuple(bam.baseName.toString(), bam) }
            .view { x -> "Pseudobulk BAM files: ${x}" }

        pseudobulk_fastqs = PSEUDOBULK_BAM_TO_FASTQ(pseudobulk_bams)

        sample_fastqs = pseudobulk_fastqs
            .map { sample_name, r1_fastq, r2_fastq ->
                tuple(sample_name, [r1_fastq], [r2_fastq])
            }
            .view { row -> "Pseudobulk sample fastqs: ${row}" }
    }

    sample_counts = NTSM_COUNT(sample_fastqs)
    NTSM_EVAL(sample_counts.collect(), modalities.collect())
}
