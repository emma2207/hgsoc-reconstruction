# Benchmarking Nextflow Pipeline

Nextflow pipeline to test tools to detect sample mislabeling (i.e. sample swaps or match samples).
This workflow is branched into two paths: one for .vcf files getting fed into tools and one for tools that need .bam files.

## Main files
- **sample-matching.yml:** yml specifying software dependencies for everything but the individual tools.
- **run_nf_pipeline.sh:** shell script to run the whole pipeline on the cluster
- **nextflow.config:** specifies global parameters and profiles for different execution environments. Note that the params can easily be overwritten via a command line argument at the time of execution.
- **main.nf:** main nextflow file. This is the file that needs to be executed to run the whole workflow. This particular main.nf is the most basic one, starting from .bam files.
    - **main_from_genotype.nf:** a variation on `main.nf` where the input is .vcf instead of .bam. I.e. the genotyping has already been performed.
    - **main_missing_samples.nf:** a variation on `main_from_genotype.nf` for pseudobulk experiments where the input data is manipulated to have missing samples.
    - **main_double_samples.nf:** a variation on `main_from_genotype.nf` for pseudobulk experiments where the input data is manipulated to have samples that are doubly represented (not literally identical samples, but drawn from the same underlying sample).
    - **main_uneven_pseudobulk_sizes.nf:** variation on `main_from_genotype.nf` for pseudobulk experiments where pseudobulk drawn from the size samples but with different numbers of cells are compared.
    - **main_pseudobulk_vs_sc.nf:** variation on `main_from_genotype.nf` comparing pseudobulks with single-cell (or single-nucleus) samples.

## Modules
Each module specifies a nextflow process that is called in main:
- **0_prepare_inputs.nf:** prepare input files .bam files
- **1a_genotype.nf:** performs variant calling for each sample and applies some (optional) filters. It results in a .vcf file and these get fed into any module whose file starts with 2a.
- **1b_filter_bam.nf:** filters the aligned reads by quality score and these filtered .bam's get fed into any module whose file starts with 2b.
- **2a_fingerprints.nf:** applies the tool CrosscheckFingerprints and needs the conda env specified in `fingerprints.yml`.
- **2a_hysys.nf:** applies the tool HaveYouSwappedYourSamples and needs the conda env specified in `hysys.yml`.
- **2a_ngscheckmate.nf:** applies the tool NGSCheckmate and needs the conda env specified in `ngscheckmate.yml`.
- **2a_vireo.nf:** applies the tool Vireo through a python script `vireo.py`.
- **2a_omicsprint.nf:** applies the Bioconductor package omicsPrint in R on SNPs extracted from the merged VCF.
- **2b_bamixchecker.mf:** applies the tool BAMixChecker needs the conda env specified in `bamixchecker.yml`.

## Parameters
Some notes on the parameter that need to be set in `nextflow.config`.

### Pseudobulk experiments
- **params.dataset:** name of dataset. Options: `hgsoc`, `high_grade_glioma`, `low_grade_glioma`, `wilms_tumor`.
- **params.pseudobulk:** is this a pseudobulk experiment or not. Must be `true` for pseudobulk experiments.
- **params.ncells:** number of cells included in the pseudobulks. Options: 10,0000; 50,000; 100,000; 500,000. 
- **params.read_depth:** the read depth filter cut-off that is used. Each dataset can define modality-specific lower thresholds plus a dataset-level `upper` cap. Options: Any number, but we focused on [0, 1, 10, 20, 30, 40].
- **params.omicsprint_max_snps:** optional cap on number of SNPs passed to omicsPrint after filtering. Use `0` to keep all SNPs.
- **params.omicsprint_call_rate:** SNP call-rate threshold passed to `omicsPrint::alleleSharing` (default in this pipeline: `0.80`; omicsPrint default: `0.95`). Lower values keep more SNPs in sparse datasets.
- **params.omicsprint_coverage_rate:** sample coverage threshold passed to `omicsPrint::alleleSharing` (default in this pipeline: `0.25`; omicsPrint default: `2/3`). Lower values retain more samples in low-coverage comparisons.

### Real data experiments
- **params.dataset:** name of dataset. Options: `hgsoc`, `hgsoc-new`, `high_grade_glioma`, `low_grade_glioma`, `wilms_tumor`.
- **params.pseudobulk:** is this a pseudobulk experiment or not. Must be `false` for real data experiments.
- **params.ncells:** number of cells included in the pseudobulks. This parameter is not used when `params.pseudobulk = false`.
- **params.read_depth:** the read depth filter cut-off that is used. Each dataset can define modality-specific lower thresholds plus a dataset-level `upper` cap. Options: Any number, but we focused on [0, 1, 10, 20, 30, 40].
- **params.omicsprint_max_snps:** optional cap on number of SNPs passed to omicsPrint after filtering. Use `0` to keep all SNPs.
- **params.omicsprint_call_rate:** SNP call-rate threshold passed to `omicsPrint::alleleSharing` (default in this pipeline: `0.80`; omicsPrint default: `0.95`). Lower values keep more SNPs in sparse datasets.
- **params.omicsprint_coverage_rate:** sample coverage threshold passed to `omicsPrint::alleleSharing` (default in this pipeline: `0.25`; omicsPrint default: `2/3`). Lower values retain more samples in low-coverage comparisons.

Additionally, for real data experiments, one has to provide the two modalities that one wants to compare. 
This is set in whichever main Nextflow files one is using (look for `modalities`).
The options are specific to each dataset:
- **hgsoc:** bulk_chunk_ribo, bulk_dissociated_polyA, bulk_dissociated_ribo & single-cell.
- **hgsoc-new:** bulk_chunk_ribo & bulk_diss_polyA.
- **high_grade_glioma:** bulk & single-cell.
- **low_grade_glioma:** bulk & single-cell.
- **wilms_tumor:** bulk & single-nucleus.
