#!/usr/bin/env nextflow

/*
 * Pipeline for sample matching analysis
 * This pipeline analyzes BAM files to determine if samples come from the same donor
 * using multiple tools: BAMixChecker, HYSYS, Vireo, NGSCheckMate, and CrossCheckFingerprints
 */

nextflow.enable.dsl = 2

// Import processes
include { PREPARE_INPUTS } from './modules/0_prepare_inputs'
include { GENOTYPE_INDIVIDUAL } from './modules/1a_genotype_individual'
include { MERGE_AND_FILTER_VCFS } from './modules/1a_merge_vcfs'
include { FILTER_BAM } from './modules/1b_filter_bam'
include { CROSSCHECK_FINGERPRINTS } from './modules/2a_fingerprints'
include { HYSYS } from './modules/2a_hysys'
include { NGSCHECKMATE } from './modules/2a_ngscheckmate'
include { VIREO_MATCH } from './modules/2a_vireo'
include { BAMIXCHECKER } from './modules/2b_bamixchecker'


// Main workflow
workflow {
    // Find all bams with the listed datatypes
    if (!params.pseudobulk) {
        modalities = ["bulk", "single-nucleus"]
        bam_patterns = modalities.collect { datatype ->
            "${params.bamsDir}/${params.dataset}/${datatype}/*/Aligned.sortedByCoord.out.bam" 
        }
        
        bam_files = channel
            .fromPath(bam_patterns)
            .map { f ->
                // safe: get the last element of the parent Path as String
                def parentName = f.parent ? f.parent.getFileName().toString() : ''
                // datatype is in the path: .../${dataset}/${datatype}/${sample}/file.bam
                def datatype = f.toString().contains("/${modalities[0]}/") ? modalities[0] :
                               f.toString().contains("/${modalities[1]}/") ? modalities[1] : 'unknown'
                tuple(parentName, f.baseName, f, datatype)
            }
            .view { x -> "Found BAM files: ${x}" }
    } else {
        modalities = ["pseudobulk"]
        bam_files = channel.fromPath("${params.pseudobulkDir}/${params.dataset}/*/pseudobulk_*_ncells_1000_*.bam")
            .map { f ->
                // safe: get the last element of the parent Path as String
                def parentName = f.parent ? f.parent.getFileName().toString() : ''
                tuple(parentName, f.baseName, f, 'pseudobulk')
            }
            .view { x -> "Found Pseudobulk BAM files: ${x}" }
    }
    
    // 0. Prepare inputs
    prepped_bam_files = PREPARE_INPUTS(bam_files)
    prepped_bam_files.view { x -> "Prepared BAM files: ${x}" }

    // 1a. Genotype individual samples in parallel
    individual_vcfs = GENOTYPE_INDIVIDUAL(prepped_bam_files)

    // 1a. Merge and filter variants
    mod_channel = channel.fromList(modalities)
    filtered_vcfs = MERGE_AND_FILTER_VCFS(
        individual_vcfs.vcf.collect(), 
        individual_vcfs.index.collect(), 
        mod_channel.collect(),
        params.pseudobulk
    )
    filtered_vcfs.modality_vcfs.view { x -> "Modality VCFs: ${x}"}
    filtered_vcfs.individual_vcfs.view{ x -> "Individual VCFs: ${x}"}

    // 1b. Filter aligned reads
    filtered_bams = FILTER_BAM(bam_files)
    filtered_bams.bam.collect().view { x -> "Filtered BAMs: ${x}" }

    // 2a. Run similarity analysis tools in parallel on filtered VCFs
    VIREO_MATCH(filtered_vcfs.modality_vcfs.collect(), mod_channel.collect())
    NGSCHECKMATE(filtered_vcfs.individual_vcfs, mod_channel.collect())
    CROSSCHECK_FINGERPRINTS(filtered_vcfs.individual_vcfs, mod_channel.collect())
    HYSYS(filtered_vcfs.individual_vcfs, mod_channel.collect())

    // 2b. Run similarity analysis tools in parallel on filtered BAMs
    BAMIXCHECKER(filtered_bams.bam.collect())
}