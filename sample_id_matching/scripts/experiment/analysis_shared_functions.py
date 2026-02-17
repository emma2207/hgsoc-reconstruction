import os
import numpy as np
import matplotlib.pyplot as plt


def all_samples_matrix_viz(
    matrix_df, pseudobulk, tool, dataset, ncells, rd, save_fig=False, FIGURES_PATH=""
):
    """
    Docstring for all_samples_matrix_viz

    :param matrix_df: Dataframe containing a heatmap / matrix of similarity measures of all samples against all samples
    :param pseudobulk: True / False
    :param tool: BAMixChecker / CrosscheckFingerprints / HYSYS / NGSCheckmate / Vireo
    :param dataset: dataset name
    :param ncells: number of cells in pseudobulk, "null" in real data
    :param rd: read depth filter cut-off
    :param sc: single-cell / single-nucleus
    :param save_fig: True / False
    :param FIGURES_PATH: path that figures get saved to
    """

    if pseudobulk:
        fig_name = f"pseudobulk_{dataset}_{tool}_rd{rd}_ncells{ncells}_all_samples_similarity_matrix"
    else:
        fig_name = f"{dataset}_{tool}_rd{rd}_all_samples_similarity_matrix"

    # Visualize ALL samples against ALL samples
    extreme_point = max(abs(matrix_df.min().min()), abs(matrix_df.max().max()))

    fig, ax = plt.subplots(figsize=(10, 10))
    if tool == "CrosscheckFingerprints":
        cmap = plt.get_cmap("RdBu")
        cmap.set_bad(color="lightgrey")
        cax = ax.matshow(matrix_df, cmap=cmap, vmin=-extreme_point, vmax=extreme_point)
    else:
        cmap = plt.get_cmap("Oranges")
        cmap.set_bad(color="lightgrey")
        cax = ax.matshow(matrix_df, cmap=cmap)

    ax.set_xticks(np.arange(len(matrix_df.columns)))
    ax.set_yticks(np.arange(len(matrix_df.index)))
    ax.set_xticklabels(
        matrix_df.columns, rotation=45, ha="left", fontdict={"fontsize": 10}
    )
    ax.set_yticklabels(matrix_df.index, fontdict={"fontsize": 10})
    ax.xaxis.set_label_position("top")
    ax.set_title(f"{tool} Similarity Matrix", pad=20)

    fig.colorbar(cax, fraction=0.046, pad=0.04, shrink=0.5)

    if save_fig:
        fig.savefig(
            os.path.join(
                FIGURES_PATH,
                f"{fig_name}.png",
            ),
            bbox_inches="tight",
            dpi=300,
        )
        fig.savefig(
            os.path.join(
                FIGURES_PATH,
                f"{fig_name}.pdf",
            ),
            bbox_inches="tight",
            dpi=300,
        )

    return


def bulk_vs_singlecell_matrix_viz(
    matrix_df,
    pseudobulk,
    tool,
    dataset,
    ncells,
    rd,
    sc="single-cell",
    save_fig=False,
    FIGURES_PATH="",
):
    """
    Docstring for bulk_vs_singlecell_matrix_viz

    :param matrix_df: Dataframe containing a heatmap / matrix of similarity measures of all samples against all samples
    :param pseudobulk: True / False
    :param tool: BAMixChecker / CrosscheckFingerprints / HYSYS / NGSCheckmate / Vireo
    :param dataset: dataset name
    :param ncells: number of cells in pseudobulk, "null" in real data
    :param rd: read depth filter cut-off
    :param sc: single-cell / single-nucleus
    :param save_fig: True / False
    :param FIGURES_PATH: path that figures get saved to
    """

    if pseudobulk:
        fig_name = f"pseudobulk_{dataset}_{tool}_rd{rd}_ncells{ncells}_bulk_vs_{sc}_similarity_matrix"
    else:
        fig_name = f"{dataset}_{tool}_rd{rd}_bulk_vs_{sc}_similarity_matrix"

    # Visualize bulk against single-cell or single-nucleus samples
    sub_matrix = matrix_df[
        [
            col
            for col in matrix_df.columns
            if "single-cell" in col or "single-nucleus" in col
        ]
    ]
    sub_matrix = sub_matrix.loc[[idx for idx in sub_matrix.index if "bulk" in idx]]

    extreme_point = max(abs(sub_matrix.min().min()), abs(sub_matrix.max().max()))

    fig, ax = plt.subplots(figsize=(8, 8))
    if tool == "CrosscheckFingerprints":
        cmap = plt.get_cmap("RdBu")
        cmap.set_bad(color="lightgrey")
        cax = ax.matshow(sub_matrix, cmap=cmap, vmin=-extreme_point, vmax=extreme_point)
    else:
        cmap = plt.get_cmap("Oranges")
        cmap.set_bad(color="lightgrey")
        cax = ax.matshow(sub_matrix, cmap=cmap)

    ax.set_xticks(np.arange(len(sub_matrix.columns)))
    ax.set_yticks(np.arange(len(sub_matrix.index)))
    ax.set_xticklabels(sub_matrix.columns, rotation=45, ha="left")
    ax.set_yticklabels(sub_matrix.index)
    ax.xaxis.set_label_position("top")
    ax.set_xlabel(f"{sc} samples")
    ax.set_ylabel("Bulk samples")
    ax.set_title(f"{tool} Bulk vs {sc} Similarity Matrix", pad=20)

    fig.colorbar(cax, fraction=0.046, pad=0.04, shrink=0.5)

    if save_fig:
        fig.savefig(
            os.path.join(
                FIGURES_PATH,
                f"{fig_name}.png",
            ),
            bbox_inches="tight",
            dpi=300,
        )
        fig.savefig(
            os.path.join(
                FIGURES_PATH,
                f"{fig_name}.pdf",
            ),
            bbox_inches="tight",
            dpi=300,
        )

    return
