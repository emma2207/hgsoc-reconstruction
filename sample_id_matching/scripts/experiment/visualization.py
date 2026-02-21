import os
import numpy as np
import matplotlib.pyplot as plt

from analysis_shared_functions import (
    parse_heatmap_matrix_bamixchecker,
    parse_heatmap_matrix_crosscheckfingerprints,
    parse_heatmap_matrix_hysys,
    parse_heatmap_matrix_ngscheckmate,
    parse_heatmap_matrix_vireo
)


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
    plt.close()

    return


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
    plt.close()

    return


def loop_heatmap_plots_pseudobulks(DATA_PATH, FIGURES_PATH, tools, datasets, ncells_list, read_depths):

    for tool in tools:
        for dataset in datasets:
            for ncells in ncells_list:
                for rd in read_depths:
                    try:
                        # Create heatmap matrix
                        if tool == "BAMixChecker":
                            _, matrix = parse_heatmap_matrix_bamixchecker(DATA_PATH, True, dataset, ncells)
                        elif tool == "CrosscheckFingerprints":
                            _, matrix = parse_heatmap_matrix_crosscheckfingerprints(DATA_PATH, True, dataset, ncells, rd)
                        elif tool == "HYSYS":
                            _, matrix = parse_heatmap_matrix_hysys(DATA_PATH, True, dataset, ncells, rd)
                        elif tool == "NGSCheckmate":
                            _, matrix = parse_heatmap_matrix_ngscheckmate(DATA_PATH, True, dataset, ncells, rd)
                        elif tool == "Vireo":
                            matrix = parse_heatmap_matrix_vireo(DATA_PATH, True, dataset, ncells, rd)
                        else:
                            print(f"Error! Do not recognize tool {tool}")
                            continue
                    except FileNotFoundError:
                        print(f"Could not find data for {tool}, pseudobulk dataset {dataset}, ncells {ncells}, read depth {rd}. ")
                        continue

                    pseudobulk_1 = "_1"
                    pseudobulk_2 = "_2"
                    matrix_filtered = matrix.loc[
                        [idx for idx in matrix.index if idx.endswith(pseudobulk_1)],
                        [col for col in matrix.columns if col.endswith(pseudobulk_2)]
                    ]

                    # Create plots
                    all_samples_matrix_viz(
                        matrix_filtered, 
                        pseudobulk=True, 
                        tool=tool, 
                        dataset=dataset, 
                        ncells=ncells, 
                        rd=rd, 
                        save_fig=True,
                        FIGURES_PATH=FIGURES_PATH
                    )


def heatmap_plot_accuracy_metrics_pseudobulk(df, dataset, metric, save_fig=False, FIGURES_PATH=""):

    fig_name = f"pseudobulk_{dataset}_{metric}_heatmap"

    fig, axes = plt.subplots(2, 2, figsize=(10, 10))
    fig.subplots_adjust(hspace=0.2)
    fig.suptitle(f'{metric} Heatmaps - {dataset} pseudobulks', fontsize=14, y=0.94)
    cmap = plt.get_cmap("Blues")
    cmap.set_bad(color="lightgrey")

    tools = ["CrosscheckFingerprints", "HYSYS", "NGSCheckmate", "Vireo"]
    for i, tool in enumerate(tools):
        sub_matrix = df[(df["dataset"] == dataset) & (df["tool"] == tool)].pivot_table(
            index=["read depth"],
            columns="ncells",
            values=metric
        )

        ax = axes[i // 2, i % 2]
        
        cax = ax.imshow(sub_matrix, cmap=cmap, vmin=0, vmax=1, aspect=1)
        
        # Add text annotations with F1 values
        for row in range(len(sub_matrix.index)):
            for col in range(len(sub_matrix.columns)):
                value = sub_matrix.iloc[row, col]
                if not np.isnan(value):
                    text_color = 'white' if value > 0.5 else 'black'
                    ax.text(col, row, f'{value:.2f}', ha='center', va='center', 
                           color=text_color, fontsize=10)
        
        ax.set_xticks(np.arange(len(sub_matrix.columns)))
        ax.set_yticks(np.arange(len(sub_matrix.index)))
        
        # Format x-tick labels in scientific notation
        xticklabels = []
        for x in sub_matrix.columns:
            exponent = int(np.log10(x))
            mantissa = x / (10 ** exponent)
            if mantissa == 1.0:
                xticklabels.append(f"$10^{{{exponent}}}$")
            else:
                xticklabels.append(f"${int(mantissa)} \\times 10^{{{exponent}}}$")
        
        ax.set_xticklabels(xticklabels, ha="center")
        ax.set_yticklabels(sub_matrix.index)
        
        # Only show x-label on bottom row
        if i // 2 == 1:
            ax.set_xlabel("# of cells")
        
        # Only show y-label on left column
        if i % 2 == 0:
            ax.set_ylabel("Read depth cut-off")
        
        ax.set_title(f"{tool}", pad=10)

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
    # fig.colorbar(cax, ax=axes.ravel().tolist(), fraction=0.046, pad=0.05, shrink=0.5)
    return