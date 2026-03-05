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
    modalities = ["pseudobulk"]
    mod_channel = channel.fromList(modalities)
    
    // Use all pseudobulk files ending in _1.vcf.gz and _2.vcf.gz
    vcf_files_1 = channel.fromPath(
    "${params.vcfsDir}/${params.dataset}/pseudobulk/ncells_${params.ncells}/*_1.vcf.gz"
    ).collect()
    vcf_files_2_all = channel.fromPath(
    "${params.vcfsDir}/${params.dataset}/pseudobulk/ncells_${params.ncells}/*_2.vcf.gz"
    )
    .collect()
    .view { n -> "Number of VCF Files 2 before removing samples: ${n.size()}"}
    index_files = channel.fromPath(
    "${params.vcfsDir}/${params.dataset}/pseudobulk/ncells_${params.ncells}/*.vcf.gz.csi"
    )

    // Remove random samples from vcf_files_2 to simulate missing data
    vcf_files_2 = vcf_files_2_all.map { files_2 ->
        def seededRandom = new Random(params.repeat) // Seed inside closure for proper serialization
        def files_2_copy = files_2.collect()
        Collections.shuffle(files_2_copy, seededRandom)
        files_2_copy.take(files_2_copy.size() - params.n_samples_to_remove)
    }
    .view { x -> "Number of VCF Files 2 after removing samples: ${x.size()}"}

    // Combine all VCF files
    vcf_files = vcf_files_1.mix(vcf_files_2)
        .collect()
        .view { files -> "All VCF files (${files.size()} total): ${files.collect { it.name }.join(', ')}" }
    
    // 1a. Filter variant calls
    filtered_vcfs = MERGE_AND_FILTER_VCFS(
        vcf_files,
        index_files.collect(),
        mod_channel.collect()
    )
    filtered_vcfs.modality_vcfs.view { x -> "VCFs by modality: ${x}" }
    filtered_vcfs.individual_vcfs.view{ x -> "Individual VCFs: ${x}"}

    // 2a. Run similarity analysis tools in parallel on filtered VCFs
    VIREO_MATCH(filtered_vcfs.modality_vcfs)
    NGSCHECKMATE(filtered_vcfs.individual_vcfs)
    CROSSCHECK_FINGERPRINTS(filtered_vcfs.individual_vcfs)
    HYSYS(filtered_vcfs.individual_vcfs, mod_channel.collect())
}