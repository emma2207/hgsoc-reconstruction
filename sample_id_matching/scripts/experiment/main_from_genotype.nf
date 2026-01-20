#!/usr/bin/env nextflow

/*
 * Pipeline for sample matching analysis
 * This pipeline analyzes BAM files to determine if samples come from the same donor
 * using multiple tools: bcftools, Vireo, NGSCheckMate, and CrossCheckFingerprints
 */

nextflow.enable.dsl = 2

// Import processes
include { FILTER_VCF } from './modules/1a_filter_vcf'
include { FILTER_BAM } from './modules/1b_filter_bam'
include { CROSSCHECK_FINGERPRINTS } from './modules/2a_fingerprints'
include { HYSYS } from './modules/2a_hysys'
include { NGSCHECKMATE } from './modules/2a_ngscheckmate'
include { VIREO_MATCH } from './modules/2a_vireo'
include { BAMIXCHECKER } from './modules/2b_bamixchecker'
include { CONPAIR } from './modules/2b_conpair'
include { BAMMATCHER } from './modules/2b_bam-matcher'


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

    // 1a. Filter variant calls
    vcf_file = channel.fromPath(
        "${params.outdir}/vcf/${params.dataset}/real_data/filtered_variants.vcf.gz"
        )
        .view { x -> "VCF File 1: ${x}"}
    mod_channel = channel.fromList(modalities)
    filtered_vcfs = FILTER_VCF(vcf_file, mod_channel.collect())
    filtered_vcfs.modality_vcfs.view { x -> "VCFs by modality: ${x}" }
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
    CONPAIR(filtered_bams.bam.collect())
    BAMMATCHER(filtered_bams.bam.collect())
}