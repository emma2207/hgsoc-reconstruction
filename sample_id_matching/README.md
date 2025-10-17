# Pseudobulk Nextflow Pipeline

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


# Sample ID Matching

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
