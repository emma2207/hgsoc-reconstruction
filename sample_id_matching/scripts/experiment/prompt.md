This is a bioinformatic Nextflow pipeline that assesses whether pairs of samples came from the same donor. The pipeline takes in aligned reads in the form of .bam files from different samples. 
1. The aligned reads get genotyped using bcftools, which results in variant calls in .vcf files. 
2. Next, quality and read depth filters are applied to the variants, also using bcftools. 
3. The filtered variant calls, are passed into several different tools that each check how similar the samples are. The existing tools that are used are Vireo, NGSCheckMate, and CrossCheckFingerprints. 
    * To use Vireo use the match_sample_id() function from match_sample_ids.py
4. The outputs from the different tools will vary, generally there will be a matrix with scores comparing all the samples against each other. We need a way to check if the results are consistent between the different tools: Do they agree on which samples are similar and which thus are likely to come from the same donor?

## Code Standards

### Nextflow structure
- Nextflow is used for the workflow. Each process should be in a separate Nextflow file.
- Bash is used in the script parts of the processes to call on the various outside tools.
- The outputs from all processes should be published in a shared directory

### Dependencies
- Include a .yml file with software dependencies that is suitable to create a conda environment with. 

## Structure
- The main file should be called experiment.nf.
- The processes will all be places in a modules/ directory.