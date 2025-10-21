import os
import pandas as pd

pl_path = "/pl/active/cgreene-sc-hgsoc/mismatch_project_data/hgsoc_ariel"

if __name__ == "__main__":
    df = pd.read_csv("../data/metadata/SraRunTable_all_hgsoc.csv")

    df = df[["Run", "Sample Name"]]
    # # Split Sample Name into Sample ID and Data Type
    temp_df = df["Sample Name"].str.split("_", n=1, expand=True)
    temp_df.columns = ["Sample ID", "Data Type"]
    df = pd.concat([df, temp_df], axis=1)

    print("Creating directories...")

    # Make directories for each data type
    data_types = list(df["Data Type"].unique())
    # for data_type in data_types:
        # os.makedirs(f"{pl_path}/{data_type}", exist_ok=True)

    print("Moving files...")

    # Move files to corresponding directories
    for _, row in df.iterrows():
        run = row["Run"]
        data_type = row["Data Type"]
        # os.system(f"mv {pl_path}/all/{run}_*.fastq.gz {pl_path}/{data_type}/")

    # os.system(f"rm -r {pl_path}/all")

    print("Finding sample matches...")

    # Figure out sample matches
    df = df[df["Data Type"]!="pooled_single_cell"]
    pivot_df = df.pivot(index="Sample ID", columns="Data Type", values="Sample Name").reset_index()
    print(pivot_df.head(10))
    pivot_df.to_csv("../data/sample_matches/sample_matches_hgsoc.csv", index=False)
    