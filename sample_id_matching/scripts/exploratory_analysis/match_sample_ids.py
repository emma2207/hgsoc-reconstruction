import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import vireoSNP


def match_sample_id(file_path1, file_path2):
    """
    Match sample IDs from dissociated bulk to bulk samples using vireoSNP.

    Args:
        file_path1 (str): Path to the first VCF file (e.g. bulk samples).
        file_path2 (str): Path to the second VCF file (e.g. dissociated bulk samples).
    Returns:
        res (dict): Dictionary containing matched sample IDs and GT probabilities.
        df (DataFrame): DataFrame with matched sample IDs.
    """

    # match bulk to dissociated bulk
    res = vireoSNP.vcf.match_VCF_samples(
        file_path1,
        file_path2,
        GT_tag1="PL",
        GT_tag2="PL",
    )

    df = pd.DataFrame(
        {
            "diss bulk sample id": res["matched_donors2"],
            "bulk sample id": res["matched_donors1"],
        }
    )
    # Clean up "sample id" columns
    df = df.replace({r"\d{6}/": ""}, regex=True)

    return res, df


def make_heatmap(cwd, df):
    """
    Create a heatmap of the 'distance' between bulk and dissociated bulk samples.
    Args:
        cwd (str): Current working directory.
        df (DataFrame): DataFrame containing 'distance matrix'.
    Returns:
        None: Saves the heatmap as a PDF and PNG file.
    """

    sorted_df = df.sort_index()
    sorted_df = sorted_df.reindex(sorted(sorted_df.columns), axis=1)

    fig, ax = plt.subplots(figsize=(15, 15))
    ax.matshow(sorted_df, cmap="Oranges")

    # Numbers in grid
    for (i, j), z in np.ndenumerate(sorted_df):
        ax.text(j, i, "{:.2f}".format(z), ha="center", va="center")

    ax.set_xticks(np.arange(len(df)))
    ax.set_yticks(np.arange(len(df)))
    ax.set_xticklabels(sorted_df.columns, rotation=45, ha="left")
    ax.set_yticklabels(sorted_df.index)
    ax.xaxis.set_label_position("top")
    ax.set_xlabel("Dissociated bulk ref")
    ax.set_ylabel("Bulk ref")
    ax.set_title("GT Prob Delta")

    filename = cwd + "/bcftools/all/GT_delta_prob_all_vs_all"
    fig.savefig(filename + ".pdf", bbox_inches="tight")
    fig.savefig(filename + ".png", bbox_inches="tight")

    return


if __name__ == "__main__":
    # Current location is /scratch/alpine/$USER/hgsoc/pooled
    cwd = os.getcwd()
    file_path1 = cwd + "/bcftools/all/bcftools_all_samples.vcf.gz"
    file_path2 = cwd + "/bcftools/all/bcftools_all_samples.vcf.gz"
    res, pair_df = match_sample_id(file_path1, file_path2)

    df = pd.DataFrame(
        res["matched_GPb_diff"],
        columns=res["matched_donors2"],
        index=res["matched_donors1"],
    )
    make_heatmap(cwd, df)

    # Save the sample pairing and the heatmap data
    pair_df.to_csv(cwd + "/bcftools/all/sample_pairing_all_v_all.csv", index=False)
    df.to_csv(cwd + "/bcftools/all/heatmap_raw_all_v_all.csv", index=True)
