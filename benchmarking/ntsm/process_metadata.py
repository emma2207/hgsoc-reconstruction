import argparse
import os

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create per-sample FASTQ mapping table for NTSM Nextflow pipeline."
    )
    parser.add_argument("--metadata", required=True, help="Path to metadata CSV")
    parser.add_argument("--fastq-dir", required=True, help="Directory containing FASTQ files")
    parser.add_argument(
        "--read-type",
        default="paired",
        choices=["single", "paired"],
        help="Read layout",
    )
    parser.add_argument("--output", required=True, help="Output CSV path")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    metadata_df = pd.read_csv(args.metadata)
    if "Run" not in metadata_df.columns or "Sample Name" not in metadata_df.columns:
        raise ValueError("Metadata must include 'Run' and 'Sample Name' columns.")

    mapping_df = metadata_df[["Run", "Sample Name"]].drop_duplicates().copy()

    if args.read_type == "paired":
        mapping_df["R1_path"] = mapping_df["Run"].apply(
            lambda run: os.path.join(args.fastq_dir, f"{run}_1.fastq.gz")
        )
        mapping_df["R2_path"] = mapping_df["Run"].apply(
            lambda run: os.path.join(args.fastq_dir, f"{run}_2.fastq.gz")
        )
    else:
        mapping_df["R1_path"] = mapping_df["Run"].apply(
            lambda run: os.path.join(args.fastq_dir, f"{run}.fastq.gz")
        )

    mapping_df["idx"] = mapping_df.groupby("Sample Name").cumcount()
    mapping_df["run_idx"] = "run_" + mapping_df["idx"].astype("Int64").astype(str)

    if args.read_type == "paired":
        output_df = mapping_df.pivot(
            index="Sample Name", columns="run_idx", values=["R1_path", "R2_path"]
        ).reset_index()
        output_df.columns = [f"{i}_{j}" if j else f"{i}" for i, j in output_df.columns]
    else:
        output_df = mapping_df.pivot(
            index="Sample Name", columns="run_idx", values="R1_path"
        ).reset_index()
        output_df.columns = [f"R1_path_{i}" if "run" in i else i for i in output_df.columns]

    output_df = output_df.dropna(axis=1, how="all").dropna(axis=0, how="all")
    output_df.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
