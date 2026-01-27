#!/usr/bin/env nextflow

/*
 * Pipeline for sample matching analysis
 * This pipeline analyzes BAM files to determine if samples come from the same donor
 * using multiple tools: bcftools, Vireo, NGSCheckMate, and CrossCheckFingerprints
 */

nextflow.enable.dsl = 2

// Import processes
include { FILTER_VCF } from './modules/1a_filter_vcf'
include { CROSSCHECK_FINGERPRINTS } from './modules/2a_fingerprints'
include { HYSYS } from './modules/2a_hysys'
include { NGSCHECKMATE } from './modules/2a_ngscheckmate'
include { VIREO_MATCH } from './modules/2a_vireo'


// Main workflow
workflow {
    modalities = ["bulk_chunk_ribo", "bulk_diss_polyA"]

    // 1. Filter variant calls
    vcf_file = channel.fromPath(
        "${params.outdir}/vcf/${params.dataset}/real_data/filtered_variants.vcf.gz"
        )
        .view { x -> "VCF File 1: ${x}"}
    mod_channel = channel.fromList(modalities)
    filtered_vcfs = FILTER_VCF(vcf_file, mod_channel.collect())
    filtered_vcfs.modality_vcfs.view { x -> "VCFs by modality: ${x}" }
    filtered_vcfs.individual_vcfs.view{ x -> "Individual VCFs: ${x}"}

    // 2. Run similarity analysis tools in parallel on filtered VCFs
    VIREO_MATCH(filtered_vcfs.modality_vcfs.collect())
    NGSCHECKMATE(filtered_vcfs.individual_vcfs)
    CROSSCHECK_FINGERPRINTS(filtered_vcfs.individual_vcfs)
    HYSYS(filtered_vcfs.individual_vcfs, mod_channel.collect())
}