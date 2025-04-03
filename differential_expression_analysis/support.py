import os
import pandas as pd

def prepare_count_matrix(INPUT_PATH: str) -> pd.DataFrame:
    """From the folder with alignment results, this function traverses
    the folder hierarchy to find all the readspergene.out.tab files and 
    combines them into a count matrix.
    """

    df = pd.DataFrame()

    for subdir, _, files in os.walk(INPUT_PATH):
        for file in files:
            if file.endswith("out.tab"):
                single_sample_counts = pd.read_csv(
                    subdir + "/" + file, 
                    delimiter="\t", 
                    header=[0, 1, 2, 3]
                )
                # Rename columns
                single_sample_counts.columns = [
                    "Gene ID", 
                    "Counts Unstranded RNA-seq", 
                    "Counts 1st Read Strand Aligned with RNA", 
                    "Counts 2nd Read Strand Aligned with RNA"
                ]
                # Grab only the 2 columns we need
                pre_count_matrix = single_sample_counts[[
                    "Gene ID", "Counts Unstranded RNA-seq"
                ]]
                # Rename the counts column with the sample name
                pre_count_matrix = pre_count_matrix.rename(
                    columns={
                        "Counts Unstranded RNA-seq": subdir[-11:]
                    }
                )
                # Merge counts from each file
                if df.empty:
                    df = pre_count_matrix
                else:
                    df = df.merge(
                        pre_count_matrix, 
                        how="outer", 
                        on="Gene ID"
                    )

    df = df.set_index("Gene ID").T
    
    return df