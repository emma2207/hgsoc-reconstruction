import argparse
import os
import pandas as pd

import vireoSNP


def match_sample_id(file1, file2, output_dir):
    """
    Match sample IDs using vireoSNP.

    Args:
        file_path (str): Path to the VCF file with all samples.
    Returns:
        res (dict): Dictionary containing matched sample IDs and GT probabilities.
        df (DataFrame): DataFrame with matched sample IDs.
    """

    # match bulk to dissociated bulk
    res = vireoSNP.vcf.match_VCF_samples(
        file1,
        file2,
        GT_tag1="PL",
        GT_tag2="PL",
    )

    df = pd.DataFrame(
        {
            "sample_id_0": res["matched_donors2"],
            "sample_id_1": res["matched_donors1"],
        }
    )

    similarity_matrix = pd.DataFrame(
        res["matched_GPb_diff"],
        index=res["matched_donors1"],
        columns=res["matched_donors2"],
    )

    # Save similarity matrix to CSV and matches
    output_dir = f"{output_dir}"
    os.makedirs(output_dir, exist_ok=True)
    similarity_matrix.to_csv(output_dir + "/similarity_matrix.csv")
    df.to_csv(output_dir + "/matched_samples.csv", index=False)

    return df, similarity_matrix


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Match sample IDs using vireoSNP"
    )
    parser.add_argument(
        "-v1",
        "--vcf_file1",
        dest="vcf_file1",
        help="Path to the VCF file with one set of samples",
        required=True,
        type=str,
    )
    parser.add_argument(
        "-v2",
        "--vcf_file2",
        dest="vcf_file2",
        help="Path to the VCF file with another set of samples",
        required=True,
        type=str,
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        dest="output_dir",
        help="Output directory name",
        required=True,
        type=str,
    )

    args = parser.parse_args()

    match_sample_id(args.vcf_file1, args.vcf_file2, args.output_dir)
