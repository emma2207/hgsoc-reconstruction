# HGSOC Reconstruction

This project has a number of different components and each has its own README.
Here are the components listed roughly by order of operations.


## (outdated) Differential Expression Analysis

Differential expression analysis between samples originating from different sites in the HGSOC dataset. 
However, this analysis was done before the realization that there were sample mislabeling issues, so the results are not to be trusted.


## Exploratory analysis

Playing around with some sample mispairing identification tools before started the benchmarking project.


## Data download

Everything regarding data download and organization for the sample mispairing benchmarking project.


## Pseudobulking (and preprocessing)

Nextflow pipeline for the preprocessing of RNA-seq data, going from raw reads (.fastq) to aligned reads (.bam). 
Pseudobulking happens at the read level by subsampling aligned reads from single-cell or single-nucleus RNA-seq.


## Benchmarking

Nextflow pipeline for the benchmarking of sample mispairing identification tools.


## Visualization

Visualization and analysis of the benchmarking results.


## Demultiplexing pooled samples

Code for the demultiplexing of pooled single-cell RNA-seq samples with the help of bulk RNA-seq references.

