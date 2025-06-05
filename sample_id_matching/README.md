# Sample ID Matching

Using various approaches to match paired bulk and dissociated bulk samples.  

## Vireo

Using Vireo's `match_VCF_samples()` function to match samples. 

- See [documentation](https://vireosnp.readthedocs.io/en/latest/index.html) for installation details.
Use either `match_sample_ids_bulk.ipynb` or `match_sample_ids.py` and `run_match_samples_py.sh` to match the samples. 
The latter option is recommended for big / slow datasets, but the functionality is the same. 
- Heatmaps of the "distance" between bulk and dissociated bulk samples are created by the scripts and `match_VCF_samples() also directly outputs the matched samples in an array.


## NGS Checkmate

Trying an available software to do the sample matching for us.

- Download and install [NGSCheckMate](https://github.com/parklab/NGSCheckMate/).
- Run NGS Checkmate using `run_ngscheckmate.sh`. It requires a .bed file that comes with the package (or you can find / make a custom one yourself), and .vcf files for individual samples, which we created using `bcftools` in `../demultiplexing_pooled_samples/scripts/split_bulk_vcf_files.sh`
- Analyze the results with `ngscheckmate_qc.ipynb`. It creates a heatmap of the correlation between the bulk and dissociated bulk samples. 

The program did not output a clear match between the bulk and the dissociated bulk samples. It is important to note that we did not play with the settings (.bed file, .vcf file filtering, etc) much, and these results might get better after optimization.
