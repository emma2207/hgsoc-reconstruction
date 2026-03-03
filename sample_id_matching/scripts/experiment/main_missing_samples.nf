#!/usr/bin/env nextflow

/*
 * Pipeline for sample matching analysis
 * This pipeline analyzes BAM files to determine if samples come from the same donor
 * using multiple tools: bcftools, Vireo, NGSCheckMate, and CrossCheckFingerprints
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
    // Find all vcfs with the listed datatypes
    n_samples_to_remove = 1
    repeats = 5
    ch_repeats = channel.from(1..repeats)
    modalities = ["pseudobulk"]
    mod_channel = channel.fromList(modalities)
    
    // Use all pseudobulk files ending in _1.vcf.gz and _2.vcf.gz
    vcf_files_1 = channel.fromPath(
    "${params.outdir}/1a_individual_vcf/${params.dataset}/pseudobulk/ncells_${params.ncells}/*_1.vcf.gz"
    ).collect()
    vcf_files_2 = channel.fromPath(
    "${params.outdir}/1a_individual_vcf/${params.dataset}/pseudobulk/ncells_${params.ncells}/*_2.vcf.gz"
    )
    .collect()
    .view { n -> "Number of VCF Files 2 before removing a sample: ${n.size()}"}
    index_files = channel.fromPath(
    "${params.outdir}/1a_individual_vcf/${params.dataset}/pseudobulk/ncells_${params.ncells}/*.vcf.gz.csi"
    )

    // Remove a random sample from vcf_files_2 to simulate missing data
    vcf_files_2 = vcf_files_2.map { 
        files -> 
        files.shuffle()
        files.take(files.size() - n_samples_to_remove)
    }
    .view { x -> "Number of VCF Files 2 after removing a sample: ${x.size()}"}

    vcf_files = vcf_files_1.mix(vcf_files_2).collect()
    .view { x -> "Total number of VCF Files: ${x.size()}"}
    
    // 1a. Filter variant calls
    // filtered_vcfs = MERGE_AND_FILTER_VCFS(
    //     vcf_files.collect(), 
    //     index_files.collect(),
    //     mod_channel.collect()
    // )
    // filtered_vcfs.modality_vcfs.view { x -> "VCFs by modality: ${x}" }
    // filtered_vcfs.individual_vcfs.view{ x -> "Individual VCFs: ${x}"}

    // // 2a. Run similarity analysis tools in parallel on filtered VCFs
    // VIREO_MATCH(filtered_vcfs.modality_vcfs.collect())
    // NGSCHECKMATE(filtered_vcfs.individual_vcfs)
    // CROSSCHECK_FINGERPRINTS(filtered_vcfs.individual_vcfs)
    // HYSYS(filtered_vcfs.individual_vcfs, mod_channel.collect())
}