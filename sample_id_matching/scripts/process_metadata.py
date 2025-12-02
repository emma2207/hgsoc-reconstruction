import argparse
import os
import pandas as pd


if __name__ == "__main__":

    # Parse command line arguments
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-d", "--dataset", dest="dataset", help="Dataset name", required=True, type=str
    )
    parser.add_argument(
        "-t",
        "--data_type",
        dest="datatype",
        help="Data type (options: all, bulk, single-cell, single-nucleus)",
        required=True,
        type=str,
    )
    parser.add_argument(
        "-r",
        "--read_type",
        dest="read_type",
        help="Read type (options: single, paired)",
        type=str,
        default="paired"
    )
    parser.add_argument(
        "-m",
        "--metadata_path",
        dest="metadata",
        help="Metadata CSV path",
        required=True,
        type=str,
    )

    args = parser.parse_args()

    # Take in Sequence Read Archive (SRA) metadata and output a run:name mapping csv
    metadata_df = pd.read_csv(args.metadata)
    mapping_df = metadata_df[["Run", "Sample Name"]].drop_duplicates()
    # Add fastq paths to the mapping dataframe
    base_path = (
        "/pl/active/cgreene-sc-hgsoc/mismatch_project_data"
        + "/" + args.dataset
        + "/" + args.datatype
    )
    if args.read_type == "paired":
        mapping_df["R1_path"] = base_path + "/" + mapping_df["Run"] + "_1.fastq.gz"
        mapping_df["R2_path"] = base_path + "/" + mapping_df["Run"] + "_2.fastq.gz"
    else:
        mapping_df["R1_path"] = base_path + "/" + mapping_df["Run"] + ".fastq.gz"
    # Count how many runs per sample and pivot the table
    mapping_df["idx"] = mapping_df.groupby("Sample Name").cumcount()
    mapping_df["run_idx"] = "run_" + mapping_df["idx"].astype("Int64").astype(str)
    if args.read_type == "paired":
        df = mapping_df.pivot(
            index="Sample Name", columns="run_idx", values=["R1_path", "R2_path"]
        ).reset_index()

        # Flatten the multi-level columns
        df.columns = [f"{i}_{j}" if j else f"{i}" for i, j in df.columns]
        df = df.dropna(axis=1, how="all").dropna(axis=0, how="all")
    else:
        df = mapping_df.pivot(
            index="Sample Name", columns="run_idx", values="R1_path"
        ).reset_index()

        df.columns = [f"R1_path_{i}" if "run" in i else i for i in df.columns]
        df = df.dropna(axis=1, how="all").dropna(axis=0, how="all")

    # Save the run:name mapping csv
    mapping_dir = "run_name_mapping"
    os.makedirs(mapping_dir, exist_ok=True)
    df.to_csv(
        mapping_dir + f"/run_name_mapping_{args.datatype}_{args.dataset}_{args.read_type}.csv",
        index=False,
    )
