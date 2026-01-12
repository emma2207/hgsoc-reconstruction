#!/usr/bin/env nextflow

/*
 * Pipeline for sample matching analysis
 * This pipeline analyzes BAM files to determine if samples come from the same donor
 * using multiple tools: bcftools, Vireo, NGSCheckMate, and CrossCheckFingerprints
 */

nextflow.enable.dsl = 2

// Import processes
include { PREPARE_INPUTS } from './modules/0_prepare_inputs'
include { GENOTYPE_AND_FILTER } from './modules/1a_genotype'
include { FILTER_BAM } from './modules/1b_filter_bam'
include { CROSSCHECK_FINGERPRINTS } from './modules/2a_fingerprints'
include { HYSYS } from './modules/2a_hysys'
include { NGSCHECKMATE } from './modules/2a_ngscheckmate'
include { VIREO_MATCH } from './modules/2a_vireo'
include { BAMIXCHECKER } from './modules/2b_bamixchecker'

// Main workflow
workflow {
    // Find all bams with the listed datatypes
    modalities = ["bulk", "single-cell"]
    bam_patterns = modalities.collect { datatype ->
        "${params.dataDir}/${params.dataset}/${datatype}/*/Aligned.sortedByCoord.out.bam" 
    }
    bam_files = channel
        .fromPath(bam_patterns)
        .map { f ->
            // safe: get the last element of the parent Path as String
            def parentName = f.parent ? f.parent.getFileName().toString() : ''
            // datatype is in the path: .../${dataset}/${datatype}/${sample}/file.bam
            def datatype = f.toString().contains('/bulk/') ? 'bulk' : 
                          f.toString().contains('/single-cell/') ? 'single-cell' : 'unknown'
            tuple(parentName, f.baseName, f, datatype)
        }
        .view { x -> "Found BAM files: ${x}" }

    // 0. Prepare inputs
    prepped_bam_files = PREPARE_INPUTS(bam_files).collect()
    prepped_bam_files.view { x -> "Prepared BAM files: ${x}" }

    // 1a. Genotype variants and filter in one step'
    mod_channel = channel.fromList(modalities)
    filtered_vcfs = GENOTYPE_AND_FILTER(prepped_bam_files, mod_channel.collect())
    filtered_vcfs.combined_vcf.view { x -> "Filtered VCFs: ${x}" }
    filtered_vcfs.modality_vcfs.view { x -> "Modality VCFs: ${x}"}
    filtered_vcfs.individual_vcfs.view{ x -> "Individual VCFs: ${x}"}

    // 1b. Filter aligned reads
    filtered_bams = FILTER_BAM(bam_files)
    filtered_bams.bam.collect().view { x -> "Filtered BAMs: ${x}" }

    // 2a. Run similarity analysis tools in parallel on filtered VCFs
    VIREO_MATCH(filtered_vcfs.modality_vcfs.collect())
    NGSCHECKMATE(filtered_vcfs.individual_vcfs)
    CROSSCHECK_FINGERPRINTS(filtered_vcfs.individual_vcfs)
    HYSYS(filtered_vcfs.individual_vcfs, mod_channel.collect())

    // 2b. Run similarity analysis tools in parallel on filtered BAMs
    BAMIXCHECKER(filtered_bams.bam.collect())
}