# Differential Expression Analysis

We looked at differentially expressed genes in bulk RNA-seq data of HGSOC, specifically looking for differentially exressed (DE) genes between different sample locations (ovary and omentum).

However, in the end we learned that the sample annotations were not reliable, and so we cannot trust the results coming out of this analysis.

## PyDESeq2

A Python implementation of the popular R package DESeq2, see [PyDESeq2 documentation](https://pydeseq2.readthedocs.io/en/stable/index.html).

The results can be found in `dea_bulk.ipynb`, and supporting functions can be found in `support_dea.py`.

## Bayexpress

I liked the idea of using a different approach to finding DE genes that does not involve "random" cut-offs. For example, performing DEA with PyDESeq2 we're setting cut-offs for the log-fold change and the adjusted p-value. We're choosing what we consider significantly differentially expressed. 

Bayesexpress, on the other hand, ranks the genes based on the Bayes factor. I thought it was a nice idea and wanted to try it. The results can be found in `bayexpress.ipynb`. In the end, I was hesitant to trust these results over PyDESeq2 because of the lack of tests on the Bayexpress code and I found the code generally poorly written. 

To get the notebook to run, you'll need supporting functions from the [Bayexpress Github](https://github.com/Morris-Research-Group/bayexpress), specifically `bayexpress_functions.py`. You can clone the repo or download that file specifically.


# Other analysis

## PCA

I was curious to see if the samples naturally cluster by sample location, so I tried some PCA in `pca.ipynb`.

## Gene-set Enrichment Analysis (GSEA)

To test the hypothesis of whether adipocyte-related genes are overexpressed in HGSOC samples from the omentum compared to samples from the ovary, we employed GSEA.
The analysis is done in `gsea.ipynb`.
