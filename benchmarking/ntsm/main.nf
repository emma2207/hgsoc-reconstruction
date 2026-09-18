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

// Single-end runs have no real R2 files; lazily create an empty gzipped
// placeholder so NTSM_COUNT's required R2 path input is always satisfied.
def emptyR2Placeholder() {
    def placeholder = file("${workDir}/ntsm_empty_R2.fastq.gz")
    if (!placeholder.exists()) {
        placeholder.parent.mkdirs()
        new java.util.zip.GZIPOutputStream(new FileOutputStream(placeholder.toString())).close()
    }
    return placeholder
}

workflow {
    // if (!params.pseudobulk) {
    //     modalities = ["bulk_dissociated_polyA", "bulk_chunk_ribo"]

    //     // Key metadata and fastq dir channels by modality, then join so each
    //     // metadata file is always paired with its own modality's fastq dir.
    //     metadata_ch = channel.fromList(modalities)
    //         .map { mod ->
    //             tuple(mod, file("${params.petaLibrary}/metadata/sra_run_tables/SraRunTable_${mod}_${params.dataset}_${params.read_type}.csv"))
    //         }
    //     fastq_dir_ch = channel.fromList(modalities)
    //         .map { mod ->
    //             tuple(mod, file("${params.petaLibrary}/raw_reads/${params.dataset}/${mod}", type: 'dir'))
    //         }
    //     metadata_ch.view { x -> "Metadata files: ${x}"}
    //     fastq_dir_ch.view { x -> "Fastq files: ${x}" }

    //     metadata_fastq_ch = metadata_ch.join(fastq_dir_ch)

    //     mapping_csv = PROCESS_METADATA(metadata_fastq_ch)

    //     sample_fastqs = mapping_csv.splitCsv(header: true)
    //             .map { row ->
    //                 def r2_paths = params.read_type == 'single' ? [emptyR2Placeholder()] : collectRunsDynamic(row, 'R2_path')
    //                 tuple(
    //                     row.'Sample Name',
    //                     collectRunsDynamic(row, 'R1_path'),
    //                     r2_paths
    //                 )
    //             }
    //         .view { row -> "Sample fastq dirs: ${row}" }
    // } 
    // else {
    //     modalities = ["pseudobulk"]
    //     pseudobulk_bam_pattern = "${params.prepared_bam_dir}/${params.dataset}/pseudobulk/ncells_${params.ncells}/*.bam"
    //     pseudobulk_bams = channel
    //         .fromPath(pseudobulk_bam_pattern)
    //         .map { bam -> tuple(bam.baseName.toString(), bam) }
    //         .view { x -> "Pseudobulk BAM files: ${x}" }

    //     pseudobulk_fastqs = PSEUDOBULK_BAM_TO_FASTQ(pseudobulk_bams)

    //     sample_fastqs = pseudobulk_fastqs
    //         .map { sample_name, r1_fastq, r2_fastq ->
    //             tuple(sample_name, [r1_fastq], [r2_fastq])
    //         }
    //         .view { row -> "Pseudobulk sample fastqs: ${row}" }
    // }
    modalities = ["bulk_chunk_ribo", "bulk_diss_polyA"]
    fastq_r1_samples = channel
        .fromPath("${params.petaLibrary}/*/*/HGSOC-*_R1_001.fastq.gz")
        .map { r1_fastq ->
            tuple(r1_fastq.parent.parent.name, r1_fastq.parent.name.take(10), r1_fastq)
        }
        .groupTuple(by: [0, 1], sort: true)

    fastq_r2_samples = channel
        .fromPath("${params.petaLibrary}/*/*/HGSOC-*_R2_001.fastq.gz")
        .map { r2_fastq ->
            tuple(r2_fastq.parent.parent.name, r2_fastq.parent.name.take(10), r2_fastq)
        }
        .groupTuple(by: [0, 1], sort: true)

    sample_fastqs = fastq_r1_samples
        .join(fastq_r2_samples, by: [0, 1])

    sample_counts = NTSM_COUNT(sample_fastqs)
    NTSM_EVAL(sample_counts.collect(), modalities.collect())
}
