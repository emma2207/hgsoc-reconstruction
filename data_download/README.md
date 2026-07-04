# Data Download
## Data

### "Original" HGSOC dataset
The dataset that the sample mismatch project, inadvertently, started with is the HGSOC dataset in the petalibrary location `/pl/active/cgreene-sc-hgsoc/ariel_sc_HGSOC`. Some notes about the data organization are in the `README_data_org.md` and the reasoning behind it can be found in [Data Detective repo](https://github.com/greenelab/hgsoc_data_detective). 

### Other datasets
Datasets that we recruited because they have bulk and single-cell RNA-seq data for matching samples and involve cancer. The are stored in the petalibrary `/pl/active/cgreene-sc-hgsoc/mismatch_project_data/$dataset/$data_type`.
- A previous (smaller) HGSOC dataset, used in [Ariel's genome biology paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-03077-7).

Datasets associated with the single-cell pediatric cancer atlas (scPCA)
- high-grade glioma
- low-grade glioma
- wilms tumor


## Downloading
All datasets in this section are available on the Sequence Read Archive (SRA). The listed scPCA-associated datasets were all publicly available, and we had to request access to the HGSOC dataset.

Looking through the scPCA database, papers, and SRA, we compiled SRA accession lists (.txt lists of SRA accession numbers) and SRA run tables (metadata .csv files) for each dataset. We renamed those files to `SRR_Acc_List_${data_type}_${dataset}.txt` and `SraRunRable_${data_type}_${dataset}.csv` respectively.

To download the data, we used the `sra-toolkit` software (available on Alpine) in a download script `download_data.sh`. The download script requires the `$dataset` and `$data_type` as inputs to find the accession lists and download the data.

Some accessions always failed to download, and we ended up downloading them manually from SRA or the European Nucleotide Archive (ENA). 

The `check_downloads.sh` script can be run to check whether all accessions were downloaded, and it outputs any missing accessions in a .txt file.
The`download_checksum.sh` verifies if any downloads are corrupted by verifying the checksum.


## Data organization
The low-grade glioma and HGSOC datasets came in a single accession list (not separated by bulk and single-cell / single-nucleus), and for downstream uses we split them into `$data_type` folders using metadata in the SRA run tables. 
`data_org_hgsoc_ariel.py` and `data_org_low_grade_glioma.py` this for the two datasets.


### Sample matches
To run the desired experiments, we need to know the presumed sample matches (i.e. which bulk RNA-seq sample should match which single-cell RNA-seq sample) for each dataset. 
For the HGSOC dataset we were able to figure this out from the SRA run table and the paper.
For all the other datasets, we contacted the data owners listed in the scPCA and asked for the information. The format in which we received the information varied and so we are processing it individually to create similar sample-match records.

- **HGSOC:** `data_org_hgsoc_ariel.py` creates the record based on the metadata.
- **High-grade glioma:** we received a csv file that was quite close to the record we wanted already, so we just did some minimal manual manipulation.
- **Low-grade glioma:** we received a csv file that was quite close to the record we wanted already, so we just did some minimal manual manipulation.
- **Wilms tumor:** `sample_matches_wilms_tumor.py` does some cleaning of the input file `SraRunTable_sn_wilms_tumor_AJMupdated.xlsx`.

The sample-match records are labelled `sample_matches_${dataset}.csv`.
