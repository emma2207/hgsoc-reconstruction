#!/usr/bin/env nextflow

/*
 * Pipeline for sample matching analysis
 * This pipeline analyzes BAM files to determine if samples come from the same donor
 * using multiple tools: bcftools, Vireo, NGSCheckMate, and CrossCheckFingerprints
 */

nextflow.enable.dsl = 2

// Import processes
include { PREPARE_INPUTS } from './modules/prepare_inputs'
include { GENOTYPE_AND_FILTER } from './modules/genotype'
include { VIREO_MATCH } from './modules/vireo'
include { NGSCHECKMATE } from './modules/ngscheckmate'
include { CREATE_FINGERPRINT_MAP } from './modules/fingerprint_map'
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
            tuple(parentName, f.baseName, f)
        }
        .view { x -> "Found BAM files: ${x}" }

    // 0. Prepare inputs
    prepped_bam_files = PREPARE_INPUTS(bam_files).collect()
    prepped_bam_files.view { x -> "Prepared BAM files: ${x}" }

    // 1a. Genotype variants and filter in one step
    filtered_vcfs = GENOTYPE_AND_FILTER(prepped_bam_files)
    filtered_vcfs.combined_vcf.view { x -> "Filtered VCFs: ${x}" }
    filtered_vcfs.individual_vcfs.view{ x -> "Individual VCFs: ${x}"}
    unique_vcf_paths = filtered_vcfs.individual_vcfs.collect().map { files -> 
        files.collect { it.parent }.unique() 
    }.flatten().view{ x -> "Folder: ${x}" }

    // 1b. Filter aligned reads
    filtered_bams = FILTER_BAM(bam_files)
    filtered_bams.bam.collect().view { x -> "Filtered BAMs: ${x}" }

    // 2a. Run similarity analysis tools in parallel on filtered VCFs
    VIREO_MATCH(filtered_vcfs.combined_vcf)
    NGSCHECKMATE(unique_vcf_paths)

    // map_file = channel.fromPath("/projects/${USER}/software/fingerprint_maps/map_files/hg38_chr.map")
    // map = CREATE_FINGERPRINT_MAP(map_file)
    CROSSCHECK_FINGERPRINTS(filtered_vcfs.individual_vcfs)
    HYSYS(filtered_vcfs.individual_vcfs)

    // // 2b. Run similarity analysis tools in parallel on filtered BAMs
    BAMIXCHECKER(filtered_bams.bam.collect())
}