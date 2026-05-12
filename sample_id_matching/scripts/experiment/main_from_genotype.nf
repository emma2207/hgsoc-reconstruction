#!/usr/bin/env nextflow

/*
 * Pipeline for sample matching analysis
 * This pipeline is a variation on main.nf that can be used if variant calling 
 * has already been performed.
 * It analyzes VCF files to determine if samples come from the same donor
 * using multiple tools: HYSYS, Vireo, NGSCheckMate, and CrossCheckFingerprints
 */

nextflow.enable.dsl = 2

// Import processes
include { MERGE_AND_FILTER_VCFS } from './modules/1a_merge_vcfs'
include { CROSSCHECK_FINGERPRINTS } from './modules/2a_fingerprints'
include { HYSYS } from './modules/2a_hysys'
include { NGSCHECKMATE } from './modules/2a_ngscheckmate'
include { VIREO_MATCH } from './modules/2a_vireo'


// Main workflow
workflow {
    // Find all bams with the listed datatypes
    if (!params.pseudobulk) {
        modalities = ["bulk", "single-cell"]
        vcf_files = channel.fromPath(
        "${params.vcfsDir}/${params.dataset}/real_data/ncells_null/*.vcf.gz"
        )
        .view { x -> "VCF Files: ${x}"}
        index_files = channel.fromPath(
        "${params.vcfsDir}/${params.dataset}/real_data/ncells_null/*.vcf.gz.csi"
        )
        .view { x -> "Index Files: ${x}"}
    } else {
        modalities = ["pseudobulk"]
        vcf_files = channel.fromPath(
        "${params.vcfsDir}/${params.dataset}/pseudobulk/ncells_${params.ncells}/*.vcf.gz"
        )
        .view { x -> "VCF Files: ${x}"}
        index_files = channel.fromPath(
        "${params.vcfsDir}/${params.dataset}/pseudobulk/ncells_${params.ncells}/*.vcf.gz.csi"
        )
        .view { x -> "Index Files: ${x}"}
    }
    mod_channel = channel.fromList(modalities)
    
    // 1a. Filter variant calls
    filtered_vcfs = MERGE_AND_FILTER_VCFS(
        vcf_files.collect(), 
        index_files.collect(),
        mod_channel.collect(),
        params.pseudobulk
    )
    filtered_vcfs.modality_vcfs.view { x -> "VCFs by modality: ${x}" }
    filtered_vcfs.individual_vcfs.view{ x -> "Individual VCFs: ${x}"}

    // 2a. Run similarity analysis tools in parallel on filtered VCFs
    VIREO_MATCH(filtered_vcfs.modality_vcfs.collect(), mod_channel.collect())
    NGSCHECKMATE(filtered_vcfs.individual_vcfs.collect(), mod_channel.collect())
    CROSSCHECK_FINGERPRINTS(filtered_vcfs.individual_vcfs.collect(), mod_channel.collect())
    HYSYS(filtered_vcfs.individual_vcfs.collect(), mod_channel.collect())
}