# Exploratory work on sample matching

Before kicking off the development of the Nextflow pipeline to test different sample mispairing detection tools, some exploratory work was done. Here are some old scripts trying out Vireo and NGSCheckMate. Note that these scripts are old and not maintained.

## Vireo
- `run_match_samples_py.sh` submits the `match_sample_ids.py` script to the HPC.
- `match_sample_ids.py` applies Vireo's `match_VCF_samples()`.
- `match_sample_ids_results.ipynb` is a notebook to analyze and visualize Vireo's results.

## NGSCheckMate
Figuring out the syntax:
 ```
 python ncm.py -V \
-bed "~/Documents/Software/NGSCheckMate/SNP/SNP_GRCh38_hg38_wChr.bed" \
–d "~/Documents/HGSOC/hgsoc-power/demultiplexing_pooled_samples/data/bcftools/sample_level_filtered/" \
–O "~/Documents/HGOSC/hgsoc-power/demultiplexing_pooled_samples/data/ngs_checkmate/" 
```
- `ngscheckmate_results.ipynb` is a notebook to analyze and visualize NGSCheckMate's results.
