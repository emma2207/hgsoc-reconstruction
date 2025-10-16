import os
import pandas as pd

pl_path = "/pl/active/cgreene-sc-hgsoc/mismatch_project_data/low_grade_glioma"

if __name__ == "__main__":
    df = pd.read_csv("../data/metadata/SraRunTable_all_low_grade_glioma.csv")

    df = df[["Run", "LibrarySource"]]

    # Move files to corresponding directories based on the data type
    for _, row in df.iterrows():
        run = row["Run"]
        library_source = row["LibrarySource"]
        data_type = "bulk" if library_source == "TRANSCRIPTOMIC" else "single-cell"
        os.system(f"mv {pl_path}/all/{run}_*.fastq.gz {pl_path}/{data_type}/")

    os.system(f"rm -r {pl_path}/all")
