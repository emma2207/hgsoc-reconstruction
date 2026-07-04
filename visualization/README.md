# Analysis of results

Analysis notebooks and supporting code to analyze the results coming out of the experimental nextflow pipeline described above.

## Analysis notebooks
- **pseudobulk_results.ipynb:** analyzes the results from the basic pseudobulk experiments. This focuses on the accuracy of the tools that we're comparing.
- **pathological_pseudobulk_experiments.ipynb:** analyzes the results from the pathological pseudobulk experiments (missing samples, double samples, uneven pseudobulk sizes, pseudobulk vs single-cell, and HYSYS heatmap + Vireo algorithm) to sanity check our results and test the limits of the tools. 
- **real_data_results.ipynb:** analyzes results from real data experiments. In any given experiment two modalities are compared and we want to check if the sample ids from those two modalities point to samples from the same patient. If not, this indicates there was a sample swap or mislabelling. 

## Supporting code
- **analysis_shared_functions.py:** has mostly data loading and manipulation stuff functions.
- **accuracy_functions.py:** has mostly functions to help calculate various accuracy metrics.
- **visualization.py:** has functions that create figures.