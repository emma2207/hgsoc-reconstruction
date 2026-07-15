# Pseudobulk Nextflow Pipeline

Nextflow pipeline to preprocess raw reads (.fastq files) and make pseudobulks *at the read level*

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