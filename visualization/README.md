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
- **standardize_heatmap_data.py:** exports selected parsed heatmap matrices as labeled, sorted CSV files.
- **standardize_sample_predictions.py:** exports labeled prediction matrices with matches as 1, non-matches as 0, and unknown results as NaN.

Use `standardize_heatmap_data.py --help` for the available options. Exports are organized by dataset, modality pair (or pseudobulk cell count), and read depth when the selected tool uses read-depth-filtered results.
By default, standardized heatmaps are stored in `../benchmarking/data/standardized_heatmaps` and predictions in `../benchmarking/data/standardized_predictions` when run from `visualization/`. Use `standardize_sample_predictions.py --help` for prediction export options. Run `standardize_heatmap_data.py` first: threshold-based tools and NTSM read those standardized heatmap CSVs, while CrosscheckFingerprints, HYSYS, NGSCheckmate, and Vireo use their explicit result files. Set `--heatmap-dir` if the standardized heatmaps are stored elsewhere. Prediction CSVs use the same directory organization. Override thresholds with repeatable `--threshold TOOL=VALUE` arguments, for example `--threshold Conpair=0.7 --threshold Peddy=0.02`; the selected threshold is recorded in the output filename.

A conda environment for everything is found in **analysis.yml**.