# Differential Expression Analysis

We looked at differentially expressed genes in bulk RNA-seq data of HGSOC, specifically looking for differentially exressed (DE) genes between different sample locations (ovary and omentum).

However, in the end we learned that the sample annotations were not reliable, and so we cannot trust the results coming out of this analysis.

In a few places, we used a list of adipocyte genetic markers from from [Emont et al, 2022](https://www.nature.com/articles/s41586-022-04518-2]), supplementary table 1.

## PyDESeq2

A Python implementation of the popular R package DESeq2, see [PyDESeq2 documentation](https://pydeseq2.readthedocs.io/en/stable/index.html).

The results can be found in `dea_bulk.ipynb`, and supporting functions can be found in `support_dea.py`.

As input to the `dea_bulk.ipynb` notebook, we mainly use a count matrix. To make the genes human-readable, we also import a mapping of gene ensemble id to gene names, and a list of adipocyte markers (Emont et al. 2022).


## Bayexpress

I liked the idea of using a different approach to finding DE genes that does not involve "random" cut-offs. 
For example, performing DEA with PyDESeq2 we're setting cut-offs for the log-fold change and the adjusted p-value. 
We're choosing what we consider significantly differentially expressed. 

Bayesexpress, on the other hand, ranks the genes based on the Bayes factor. I thought it was a nice idea and wanted to try it. The results can be found in `bayexpress.ipynb`. 
In the end, I was hesitant to trust these results over PyDESeq2 because of the lack of tests on the Bayexpress code and I found the code generally poorly written. 

To get the notebook to run, you'll need supporting functions from the [Bayexpress Github](https://github.com/Morris-Research-Group/bayexpress), specifically `bayexpress_functions.py`. 
You can clone the repo or download that file specifically.

Similar the the pydeseq2 analysis, the `bayexpress.ipynb` notebook, needs a count matrix, and a mapping of gene ensemble id to gene names. 
Additionally, for comparison with the pyDESeq2 analysis, we use a list of differentially expressed genes as found by the pyDESeq2 analysis.


# Other analysis

## PCA

I was curious to see if the samples naturally cluster by sample location, so I tried some PCA in `pca.ipynb`.

## Gene-set Enrichment Analysis (GSEA)

To test the hypothesis of whether adipocyte-related genes are overexpressed in HGSOC samples from the omentum compared to samples from the ovary, we employed GSEA. We used the [gseapy package](https://gseapy.readthedocs.io/en/latest/index.html). 

The analysis is done in `gsea.ipynb` and uses as inputs a list of differentially expressed genes (output from `dea_bulk.ipynb`) and a list of adipocyte genetic markers (Emont et al. 2022).

The main output is a plot of the results showing how the adipocyte markers are distributed throughout the ranked genes.
The genes are ranked by a combination of the adjusted p-value and the log-fold change. 
