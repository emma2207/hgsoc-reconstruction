#!/usr/bin/env nextflow

/*
 * Pipeline for sample matching analysis
 * This pipeline analyzes BAM files to determine if samples come from the same donor
 * using multiple tools: bcftools, Vireo, NGSCheckMate, and CrossCheckFingerprints
 */

nextflow.enable.dsl = 2

// Import processes
include { GENOTYPE_AND_FILTER } from './modules/genotype'
include { VIREO_MATCH } from './modules/vireo'
include { NGSCHECKMATE } from './modules/ngscheckmate'
include { CROSSCHECK_FINGERPRINTS } from './modules/fingerprints'
include { COMPARE_RESULTS } from './modules/compare'

// Main workflow
workflow {
    // Create channel from BAM files
    bam_files = channel
        .fromPath(params.bams)
        // get parent folder name from a Path
        .map { f ->
            // safe: get the last element of the parent Path as String
            def parentName = f.parent ? f.parent.getFileName().toString() : ''
            tuple(parentName, f)
        }
        .view { x -> "Found BAM files: ${x}" }

    // 1. Genotype variants and filter in one step
    filtered_vcfs = GENOTYPE_AND_FILTER(bam_files)

    // 2. Run similarity analysis tools in parallel
    vireo_results = VIREO_MATCH(filtered_vcfs)
    ngscheckmate_results = NGSCHECKMATE(filtered_vcfs)
    fingerprint_results = CROSSCHECK_FINGERPRINTS(filtered_vcfs)

    // 4. Compare results from different tools
    COMPARE_RESULTS(
        vireo_results,
        ngscheckmate_results,
        fingerprint_results
    )
}