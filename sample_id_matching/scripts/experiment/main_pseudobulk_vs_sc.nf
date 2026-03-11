#!/usr/bin/env nextflow

/*
 * Pipeline for sample matching analysis
 * This pipeline analyzes BAM files to determine if samples come from the same donor
 * using multiple tools: bcftools, Vireo, NGSCheckMate, and CrossCheckFingerprints
 */

nextflow.enable.dsl = 2

// Import processes
include { MERGE_AND_FILTER_VCFS as MERGE_PSEUDOBULK_VCFS} from './modules/1a_merge_vcfs'
include { MERGE_AND_FILTER_VCFS as MERGE_SC_VCFS} from './modules/1a_merge_vcfs'
include { CROSSCHECK_FINGERPRINTS } from './modules/2a_fingerprints'
include { HYSYS } from './modules/2a_hysys'
include { NGSCHECKMATE } from './modules/2a_ngscheckmate'
include { VIREO_MATCH } from './modules/2a_vireo'


// Main workflow
workflow {
    sc_vcf_files = channel.fromPath(
    "${params.vcfsDir}/${params.dataset}/real_data/ncells_null/*_single-cell.vcf.gz"
    ).collect()
    .view { x -> "SC VCF Files: ${x}"}
    sc_index_files = channel.fromPath(
    "${params.vcfsDir}/${params.dataset}/real_data/ncells_null/*_single-cell.vcf.gz.csi"
    ).collect()

    pseudobulk_vcf_files = channel.fromPath(
    "${params.vcfsDir}/${params.dataset}/pseudobulk/ncells_${params.ncells}/*_1.vcf.gz"
    ).collect()
    .view { x -> "Pseudobulk VCF Files: ${x}"}
    pseudobulk_index_files = channel.fromPath(
    "${params.vcfsDir}/${params.dataset}/pseudobulk/ncells_${params.ncells}/*_1.vcf.gz.csi"
    ).collect()
    
    // 1a. Filter variant calls
    pseudobulk_channel = channel.from(["pseudobulk"])
    pseudobulk_channel.view { x -> "Pseudobulk channel: ${x}" }
    filtered_pseudobulk_vcfs = MERGE_PSEUDOBULK_VCFS(
        pseudobulk_vcf_files, 
        pseudobulk_index_files,
        pseudobulk_channel,
    )
    filtered_pseudobulk_vcfs.individual_vcfs.view{ x -> "Individual VCFs: ${x}"}

    sc_channel = channel.from(["single-cell"])
    sc_channel.view { x -> "Single-cell channel: ${x}" }
    filtered_sc_vcfs = MERGE_SC_VCFS(
        sc_vcf_files, 
        sc_index_files,
        sc_channel,
    )
    filtered_sc_vcfs.individual_vcfs.view{ x -> "Individual VCFs: ${x}"}

    modality_vcfs = filtered_pseudobulk_vcfs.modality_vcfs.mix(filtered_sc_vcfs.modality_vcfs).collect()
    individual_vcfs = filtered_pseudobulk_vcfs.individual_vcfs.mix(filtered_sc_vcfs.individual_vcfs).collect()
    modality_vcfs.view { x -> "VCFs by modality: ${x}" }
    individual_vcfs.view{ x -> "Individual VCFs: ${x}"}

    mod_channel = pseudobulk_channel.mix(sc_channel).collect()
    // 2a. Run similarity analysis tools in parallel on filtered VCFs
    VIREO_MATCH(modality_vcfs)
    NGSCHECKMATE(individual_vcfs)
    CROSSCHECK_FINGERPRINTS(individual_vcfs)
    HYSYS(individual_vcfs, mod_channel)
}