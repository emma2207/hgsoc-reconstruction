import argparse
import os
import pandas as pd
import re


def create_output_dirs(base_dir, dataset, sample_names):
    for sample in sample_names:
        dir_path = os.path.join(base_dir, 'fastq', dataset, sample)
        os.makedirs(dir_path, exist_ok=True)
    return 


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
        "-m",
        "--metadata_path",
        dest="metadata",
        help="Metadata CSV path",
        required=True,
        type=str,
    )

    args = parser.parse_args()

    # Take in SRA metadata and output a run:name mapping csv
    metadata_df = pd.read_csv(args.metadata)
    mapping_df = metadata_df[["Run", "Sample Name"]].drop_duplicates()
    mapping_df.to_csv(
        f"../data/run-name_mapping/run_name_mapping_{args.datatype}_{args.dataset}.csv", index=False
    )

    # Use dataset and sample names to create a directory structure
    base_dir = os.getcwd() + "/../data"
    create_output_dirs(
        base_dir=base_dir,
        dataset=args.dataset,
        sample_names=mapping_df["Sample Name"].unique(),
    )

    # Use run:name mapping to move fastq.gz files to appropriate directory
    # (Organized by dataset, datatype, sample name)
    sample_dir = os.path.join(base_dir, "fastq", args.dataset, args.datatype)
    srr_files = [f for f in os.listdir(sample_dir) if f.endswith(".fastq.gz")]

    for srr in srr_files:
        run_id = re.sub(r"_\d{1}_truncated.fastq.gz", "", srr)
        sample_name = mapping_df[mapping_df["Run"] == run_id]["Sample Name"].values
        if len(sample_name) == 1:
            sample_name = sample_name[0]
            dest_path = os.path.join(base_dir, "fastq", args.dataset, args.datatype, sample_name)
            os.makedirs(dest_path, exist_ok=True)
            os.rename(
                os.path.join(sample_dir, srr),
                os.path.join(dest_path, srr),
            )
        else:
            print(f"Warning: Run ID {run_id} has {len(sample_name)} associated samples.")