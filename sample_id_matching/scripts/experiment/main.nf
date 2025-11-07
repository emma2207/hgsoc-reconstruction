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
include { HYSYS } from './modules/hysys'
include { FILTER_BAM } from './modules/filter_bam'
include { BAMIXCHECKER } from './modules/bamixchecker'

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

    // 1a. Genotype variants and filter in one step
    filtered_vcfs = GENOTYPE_AND_FILTER(bam_files)

    // 1b. Filter aligned reads
    filtered_bams = FILTER_BAM(bam_files)

    // 2a. Run similarity analysis tools in parallel on filtered VCFs
    VIREO_MATCH(filtered_vcfs)
    NGSCHECKMATE(filtered_vcfs)
    CROSSCHECK_FINGERPRINTS(filtered_vcfs)
    HYSYS(filtered_vcfs)

    // 2b. Run similarity analysis tools in parallel on filtered BAMs
    BAMIXCHECKER(filtered_bams)
}