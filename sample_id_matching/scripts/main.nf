#!/usr/bin/env nextflow

/*
 * main.nf — run pseudobulking or bulk preprocessing
 */

nextflow.enable.dsl = 2

include { PROCESS_METADATA } from './modules/process_metadata.nf'
include { QC_READS_WITH_FASTP } from './modules/fastp.nf'
include { ALIGNMENT_WITH_STAR } from './modules/star_alignment.nf'
include { CREATE_PSEUDOBULKS } from './modules/make_pseudobulks.nf'

params.debug = params.get('debug', false)

// Helper to collect R1/R2 files for any matching run columns,
// preserving numeric order of run indices.
def collectRunsDynamic(row, prefix) {
    // collect keys and build a plain-string regex (no slashy regex interpolation)
    def keys   = row.keySet().toList()
    def regex  = "${prefix}_run_\\d+"        // note: double-escape for \d inside a GString
    def matched = keys.findAll { k -> k ==~ regex }
    def sorted = matched.sort { a, b ->
        (a =~ /\d+$/)[0].toInteger() <=> (b =~ /\d+$/)[0].toInteger()
    }
    def values = sorted.collect { k -> row[k]?.toString()?.trim() }.findAll { it }
    def files = values.collect { file(it) }

    return files
}

workflow {
    // Input metadata & fastqs
    metadata = channel.fromPath("${params.petaLibrary}/metadata/sra_run_tables/SraRunTable_${params.datatype}_${params.dataset}.csv")
    input_data = channel.fromPath("${params.petaLibrary}/${params.dataset}/${params.datatype}", type: 'dir')

    // Process metadata
    fastqs = PROCESS_METADATA(metadata, input_data)

    // Split the fastqs by R1 and R2 paths for each sample
    sample_ch = fastqs.splitCsv(header: true)
            .map { row ->
                tuple(
                    row.'Sample Name',
                    collectRunsDynamic(row, 'R1_path'),
                    collectRunsDynamic(row, 'R2_path')
                )
            }
        .view { row -> "Sample fastq dirs: ${row}" }

    // QC
    fastp_out = QC_READS_WITH_FASTP(sample_ch)
    fastp_out.merged_fastqs.view { x -> "Fastp merged reads: ${x}" }
    
    // Align fastqs
    bam = ALIGNMENT_WITH_STAR(fastp_out.merged_fastqs)
    bam.aligned_reads.view { x -> "Aligned Reads: ${x}" }

    // Create pseudobulks
    if (params.datatype == 'single-cell' || params.datatype == 'single-nucleus') {
        println "Creating pseudobulks for datatype: ${params.datatype}"
        
        n_barcodes = channel.of(params.n_barcodes)
        n_pseudobulks = channel.of(params.n_pseudobulks)
        pseudobulk_input = bam.aligned_reads.combine(n_barcodes).combine(n_pseudobulks)
        CREATE_PSEUDOBULKS(pseudobulk_input)
    } else {
        println "Skipping pseudobulk creation for datatype: ${params.datatype}"
        return
    }
    
}
