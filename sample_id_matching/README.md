# 1. Data organization

## "Original" HGSOC dataset
The dataset that the sample mismatch project, inadvertently, started with is the HGSOC dataset in the petalibrary location `/pl/active/cgreene-sc-hgsoc/ariel_sc_HGSOC`. This PR does not touch that data. It only adds some notes about the data organization in the `README_data_org.md`. This README was also added to the petalibrary folder.

## Other datasets
Datasets that we recruited because they have bulk and single-cell RNA-seq data for matching samples and involve cancer. The are stored in the petalibrary `/pl/active/cgreene-sc-hgsoc/mismatch_project_data/$dataset/$data_type`.
- A previous (smaller) HGSOC dataset, used in [Ariel's genome biology paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-03077-7).

Datasets associated with the single-cell pediatric cancer atlas (scPCA)
- central nervous system tumors
- high-grade glioma
- low-grade glioma
- wilms tumor

### Data download
See `/scripts/data_download`.

All datasets in this section are available on the Sequence Read Archive (SRA). The listed scPCA-associated datasets were all publicly available, and we had to request access to the HGSOC dataset.

Looking through the scPCA database, papers, and SRA, we compiled SRA accession lists (.txt lists of SRA accession numbers) and SRA run tables (metadata .csv files) for each dataset. We renamed those files to `SRR_Acc_List_${data_type}_${dataset}.txt` and `SraRunRable_${data_type}_${dataset}.csv` respectively.

To download the data, we used the `sra-toolkit` software (available on Alpine) in a download script `download_data.sh`. The download script requires the `$dataset` and `$data_type` as inputs to find the accession lists and download the data.

Some accessions always failed to download, and we ended up downloading them manually from SRA or the European Nucleotide Archive (ENA). 

The `check_downloads.sh` script can be run to check whether all accessions were downloaded, and it outputs any missing accessions in a .txt file.

### Data organization
See `/scripts/data_download`.

The low-grade glioma and HGSOC datasets came in a single accession list (not separated by bulk and single-cell / single-nucleus), and for downstream uses we split them into `$data_type` folders using metadata in the SRA run tables. 
`data_org_hgsoc_ariel.py` and `data_org_low_grade_glioma.py` this for the two datasets.

### Sample matches
See `/scripts/data_download`.

To run the desired experiments, we need to know the presumed sample matches (i.e. which bulk RNA-seq sample should match which single-cell RNA-seq sample) for each dataset. 
For the HGSOC dataset we were able to figure this out from the SRA run table and the paper.
For all the other datasets, we contacted the data owners listed in the scPCA and asked for the information. The format in which we received the information varied and so we are processing it individually to create similar sample-match records.

- **HGSOC:** `data_org_hgsoc_ariel.py` creates the record based on the metadata.
- **Central nervous system tumors:** `sample_matches_cnst.py` required quite a bit of munging and is based on a spreadsheet `sample_matches_central_nervous_system_tumor.xlsx` received from the data owner.
- **High-grade glioma:** we received a csv file that was quite close to the record we wanted already, so we just did some minimal manual manipulation.
- **Low-grade glioma:** we received a csv file that was quite close to the record we wanted already, so we just did some minimal manual manipulation.
- **Wilms tumor:** `sample_matches_wilms_tumor.py` does some cleaning of the input file `SraRunTable_sn_wilms_tumor_AJMupdated.xlsx`.

The sample-match records are labelled `sample_matches_${dataset}.csv`.

### Pool contents
A few datasets contain pooled single-cell / single-nucleus data. We want a record of which samples are contained in each pool.

- **HGSOC:** Pulled the information from the paper.
- **Central nervous system tumors:** the information is pulled from the `sample_matches_central_nervous_system_tumor.xlsx` spreadsheet in the `sample_matches_cnst.py` script.

The pool-contents records are labelled `pool_contents_${dataset}.csv`.


# 2. Exploratory Analysis Sample ID Matching

Using various approaches to match paired bulk and dissociated bulk samples.  

## Vireo

Using Vireo's `match_VCF_samples()` function to match samples. 

- See [documentation](https://vireosnp.readthedocs.io/en/latest/index.html) for installation details.
Use either `match_sample_ids_bulk.ipynb` or `match_sample_ids.py` and `run_match_samples_py.sh` to match the samples. 
The latter option is recommended for big / slow datasets, but the functionality is the same. 
- Heatmaps of the "distance" between bulk and dissociated bulk samples are created by the scripts and `match_VCF_samples()` also directly outputs the matched samples in an array.

## NGS Checkmate

Trying an available software to do the sample matching for us.

- Download and install [NGSCheckMate](https://github.com/parklab/NGSCheckMate/).
- Run NGS Checkmate using `run_ngscheckmate.sh`. I actually ran the command directly in the terminal, but I wanted to keep an example of the syntax. It requires a .bed file that comes with the package (or you can find / make a custom one yourself), and .vcf files for individual samples, which we created using `bcftools` in `../demultiplexing_pooled_samples/scripts/split_bulk_vcf_files.sh`
- Analyze the results with `ngscheckmate_qc.ipynb`. It creates a heatmap of the correlation between the bulk and dissociated bulk samples. 

The program did not output a clear match between the bulk and the dissociated bulk samples. It is important to note that we did not play with the settings (.bed file, .vcf file filtering, etc) much, and these results might get better after optimization.


# 3. Pseudobulk Nextflow Pipeline
See `/scripts/pseudobulking`.

Nextflow pipeline to make pseudobulks *at the read level*

## Main files
- **simulation_environment.yml:** yml specifying software dependencies. Note that I have not been able to install subset-bam on my local computer, I have only been able to use it on the cluster. 
- **run_nf_pipeline.sh:** shell script to run the whole pipeline on the cluster
- **nextflow.config:** specifies global parameters and profiles for different execution environments. Note that the params can easily be overwritten via a command line argument at the time of execution.
- **main.nf:** main nextflow file. This is the file that needs to be execute to execute the whole workflow. 

## Modules

Each module specifies a nextflow process that is called in main. 
In order of operations:
- **process_metadata.nf:** processes metadata files using process_metadata.py
    - **process_metadata.py:** python script that takes in metadata taken from the SRA Run Selector, figures out which SRA runs correspond to which samples, and writes the names and paths to a csv file
- **fastp.nf:** QC of fastq files using fastp
- **star_alignment.nf:** alignes reads from fastq files for a single sample with a reference genome using STAR
- **make_pseudobulks.nf** creates pseudobulks by randomly selecting a number of cell barcodes using a python script and then creating .bam files with all the reads corresponding to those cells using subset-bam
    - **select_barcodes.py:** takes aligned reads from single-cell or single-nucleus RNA-seq, finds the unique barcodes and randomly samples them.


# 4. Experimental Nextflow Pipeline
See `/scripts/experiment`.

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
- **prompt.md:** original Github copilot prompt to get the basics of the pipeline set up. The current version of the pipeline is very different from the pipeline resulting directly from this prompt.

## Modules
Each module specifies a nextflow process that is called in main:
- **0_prepare_inputs.nf:** prepare input files .bam files
- **1a_genotype.nf:** performs variant calling for each sample and applies some (optional) filters. It results in a .vcf file and these get fed into any module whose file starts with 2a.
- **1b_filter_bam.nf:** filters the aligned reads by quality score and these filtered .bam's get fed into any module whose file starts with 2b.
- **2a_fingerprints.nf:** applies the tool CrosscheckFingerprints and needs the conda env specified in `fingerprints.yml`.
- **2a_hysys.nf:** applies the tool HaveYouSwappedYourSamples and needs the conda env specified in `hysys.yml`.
- **2a_ngscheckmate.nf:** applies the tool NGSCheckmate and needs the conda env specified in `ngscheckmate.yml`.
- **2a_vireo.nf:** applies the tool Vireo through a python script `vireo.py`.
- **2b_bamixchecker.mf:** applies the tool BAMixChecker needs the conda env specified in `bamixchecker.yml`.

## Parameters
Some notes on the parameter that need to be set in `nextflow.config`.

### Pseudobulk experiments
- **params.dataset:** name of dataset. Options: `hgsoc`, `high_grade_glioma`, `low_grade_glioma`, `wilms_tumor`.
- **params.pseudobulk:** is this a pseudobulk experiment or not. Must be `true` for pseudobulk experiments.
- **params.ncells:** number of cells included in the pseudobulks. Options: 10,0000; 50,000; 100,000; 500,000. 
- **params.read_depth:** the read depth filter cut-off that is used. Each dataset can define modality-specific lower thresholds plus a dataset-level `upper` cap. Options: Any number, but we focused on [0, 1, 10, 20, 30, 40].

### Real data experiments
- **params.dataset:** name of dataset. Options: `hgsoc`, `hgsoc-new`, `high_grade_glioma`, `low_grade_glioma`, `wilms_tumor`.
- **params.pseudobulk:** is this a pseudobulk experiment or not. Must be `false` for real data experiments.
- **params.ncells:** number of cells included in the pseudobulks. This parameter is not used when `params.pseudobulk = false`.
- **params.read_depth:** the read depth filter cut-off that is used. Each dataset can define modality-specific lower thresholds plus a dataset-level `upper` cap. Options: Any number, but we focused on [0, 1, 10, 20, 30, 40].

Additionally, for real data experiments, one has to provide the two modalities that one wants to compare. 
This is set in whichever main Nextflow files one is using (look for `modalities`).
The options are specific to each dataset:
- **hgsoc:** bulk_chunk_ribo, bulk_dissociated_polyA, bulk_dissociated_ribo & single-cell.
- **hgsoc-new:** bulk_chunk_ribo & bulk_diss_polyA.
- **high_grade_glioma:** bulk & single-cell.
- **low_grade_glioma:** bulk & single-cell.
- **wilms_tumor:** bulk & single-nucleus.


# 5. Analysis of results
See `/scripts/experiment`.

Analysis notebooks and supporting code to analyze the results coming out of the experimental nextflow pipeline described above.

## Analysis notebooks
- **pseudobulk_results.ipynb:** analyzes the results from the basic pseudobulk experiments. This focuses on the accuracy of the tools that we're comparing.
- **pathological_pseudobulk_experiments.ipynb:** analyzes the results from the pathological pseudobulk experiments (missing samples, double samples, uneven pseudobulk sizes, pseudobulk vs single-cell, and HYSYS heatmap + Vireo algorithm) to sanity check our results and test the limits of the tools. 
- **real_data_results.ipynb:** analyzes results from real data experiments. In any given experiment two modalities are compared and we want to check if the sample ids from those two modalities point to samples from the same patient. If not, this indicates there was a sample swap or mislabelling. 

## Supporting code
- **analysis_shared_functions.py:** has mostly data loading and manipulation stuff functions.
- **accuracy_functions.py:** has mostly functions to help calculate various accuracy metrics.
- **visualization.py:** has functions that create figures.