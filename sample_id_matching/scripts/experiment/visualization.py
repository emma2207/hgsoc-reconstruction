import os
import numpy as np
import matplotlib.pyplot as plt

from analysis_shared_functions import (
    load_heatmap_data,
    order_matrix_by_expected_matches,
)
from accuracy_functions import (
    load_expected_matches_real_data,
    accuracy_metrics_averaged_over_iterations,
)


def bulk_vs_singlecell_matrix_viz(
    DATA_PATH,
    pseudobulk,
    tool,
    dataset,
    ncells,
    rd,
    mod1="bulk",
    mod2="single-cell",
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
    :param mod1: first modality (e.g., "bulk")
    :param mod2: second modality (e.g., "single-cell")
    :param save_fig: True / False
    :param FIGURES_PATH: path that figures get saved to
    """

    if pseudobulk:
        fig_name = f"pseudobulk_{dataset}_{tool}_rd{rd}_ncells{ncells}_{mod1}_vs_{mod2}_similarity_matrix"
    else:
        fig_name = f"{dataset}_{tool}_rd{rd}_{mod1}_vs_{mod2}_similarity_matrix"

    # Visualize bulk against single-cell or single-nucleus samples
    matrix_df = load_heatmap_data(DATA_PATH, pseudobulk, tool, dataset, ncells, rd, mod1, mod2)    
    if matrix_df is None:
        print(f"No matrix data available for {tool} on {dataset} with read depth {rd} and ncells {ncells}. Skipping plot.")
        return
    if not pseudobulk:
        expected_matches = load_expected_matches_real_data(
            DATA_PATH, 
            dataset
        )
        matrix_df = order_matrix_by_expected_matches(
            matrix_df,
            expected_matches,
            tool,
            mod1,
            mod2,
        )
    if matrix_df is not None:
        sub_matrix = matrix_df[[col for col in matrix_df.columns if mod2 in col]]
        sub_matrix = sub_matrix.loc[[idx for idx in sub_matrix.index if mod1 in idx]]

        extreme_point = max(abs(sub_matrix.min().min()), abs(sub_matrix.max().max()))

        fig, ax = plt.subplots(figsize=(8, 8))
        if tool == "CrosscheckFingerprints":
            cmap = plt.get_cmap("RdBu")
            cmap.set_bad(color="lightgrey")
            cax = ax.matshow(
                sub_matrix, cmap=cmap, vmin=-extreme_point, vmax=extreme_point
            )
        else:
            cmap = plt.get_cmap("Oranges")
            cmap.set_bad(color="lightgrey")
            cax = ax.matshow(sub_matrix, cmap=cmap)

        ax.set_xticks(np.arange(len(sub_matrix.columns)))
        ax.set_yticks(np.arange(len(sub_matrix.index)))
        ax.set_xticklabels(sub_matrix.columns, rotation=45, ha="left")
        ax.set_yticklabels(sub_matrix.index)
        ax.xaxis.set_label_position("top")
        ax.set_xlabel(f"{mod2} samples")
        ax.set_ylabel(f"{mod1} samples")
        ax.set_title(f"{tool} {mod1} vs {mod2} similarity matrix", pad=20)

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


def all_samples_matrix_viz(
    DATA_PATH, pseudobulk, tool, dataset, ncells, rd, save_fig=False, FIGURES_PATH="", mod1="bulk", mod2="single-cell"
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

    matrix_df = load_heatmap_data(DATA_PATH, pseudobulk, tool, dataset, ncells, rd, mod1, mod2)
    if matrix_df is not None:
        # Visualize ALL samples against ALL samples
        extreme_point = max(abs(matrix_df.min().min()), abs(matrix_df.max().max()))

        fig, ax = plt.subplots(figsize=(10, 10))
        if tool == "CrosscheckFingerprints":
            cmap = plt.get_cmap("RdBu")
            cmap.set_bad(color="lightgrey")
            cax = ax.matshow(
                matrix_df, cmap=cmap, vmin=-extreme_point, vmax=extreme_point
            )
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
        ax.set_title(f"{tool} similarity matrix", pad=20)

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


def super_plot_heatmaps_pseudobulk(
    DATA_PATH,
    FIGURES_PATH,
    tool,
    dataset,
    ncells_list,
    read_depths,
    experiment,
    save_fig=False,
):

    fig_name = f"superplot_pseudobulk_{dataset}_{tool}_all_samples_similarity_matrix"
    super_extreme_point = 0
    all_matrices = []

    # Find all the data for the plots
    for ncells in ncells_list:
        for rd in read_depths:
            matrix = load_heatmap_data(DATA_PATH, True, tool, dataset, ncells, rd)

            if matrix is not None:
                if experiment == "pseudobulk_vs_sc":
                    pseudobulk_1 = "_1"
                    pseudobulk_2 = "_single-cell"
                else:
                    pseudobulk_1 = "_1"
                    pseudobulk_2 = "_2"
                matrix_filtered = matrix.loc[
                    [idx for idx in matrix.index if idx.endswith(pseudobulk_1)],
                    [col for col in matrix.columns if col.endswith(pseudobulk_2)],
                ]

                # Store matrix for later use in super plot
                all_matrices.append(
                    {"ncells": ncells, "rd": rd, "matrix": matrix_filtered}
                )

                # Calculate the extreme point over all matrices to use the same color scale for all heatmaps in the super plot
                extreme_point = max(
                    abs(matrix_filtered.min().min()), abs(matrix_filtered.max().max())
                )
                if extreme_point > super_extreme_point:
                    super_extreme_point = extreme_point

    # Loop through the data again to create the super plot
    fig, axes = plt.subplots(
        len(read_depths),
        len(ncells_list),
        figsize=(5 * len(ncells_list), 4 * len(read_depths)),
        sharex=True,
        sharey=True,
    )
    fig.subplots_adjust(hspace=0.05, wspace=0.05, top=0.95)
    if len(read_depths) == 1:
        yval = 1.1
    else:
        yval = 0.98
    fig.suptitle(
        f"{tool} - {dataset} pseudobulks similarity matrices", fontsize=16, y=yval
    )

    for ncells in ncells_list:
        for rd in read_depths:
            if len(read_depths) == 1:
                ax = axes[ncells_list.index(ncells)]
            elif len(ncells_list) == 1:
                ax = axes[read_depths.index(rd)]
            else:
                ax = axes[read_depths.index(rd), ncells_list.index(ncells)]
            matrix_list = [
                info["matrix"]
                for info in all_matrices
                if info["ncells"] == ncells and info["rd"] == rd
            ]

            if not matrix_list or matrix_list[0] is None:
                # Skip this subplot if no data available
                ax.axis("off")
                ax.text(
                    0.5,
                    0.5,
                    "No data",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                    fontsize=12,
                )
                continue

            matrix = matrix_list[0]

            if tool == "CrosscheckFingerprints":
                cmap = plt.get_cmap("RdBu")
                cmap.set_bad(color="lightgrey")
                cax = ax.imshow(
                    matrix,
                    cmap=cmap,
                    vmin=-super_extreme_point,
                    vmax=super_extreme_point,
                )
            else:
                cmap = plt.get_cmap("Oranges")
                cmap.set_bad(color="lightgrey")
                cax = ax.imshow(matrix, cmap=cmap)
            ax.set_xticks(np.arange(len(matrix.columns)))
            ax.set_yticks(np.arange(len(matrix.index)))
            ax.set_xticklabels([""] * len(matrix.columns))
            ax.set_yticklabels([""] * len(matrix.index))
            # Add title with ncells to the top row
            if rd == read_depths[0]:
                ax.set_title(f"{ncells} cells", pad=10, fontsize=12)
            # Add text with read depth to the right column
            if ncells == ncells_list[-1]:
                ax.text(
                    1.05,
                    0.5,
                    f"Read depth filter: {rd}",
                    transform=ax.transAxes,
                    fontsize=12,
                    rotation=90,
                    va="center",
                )

    # Add one common colorbar for all subplots
    # cbar = fig.colorbar(
    #     cax, ax=axes.ravel().tolist(), fraction=0.046, pad=0.05, shrink=0.5
    # )

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


def loop_heatmap_plots_real_data(
    DATA_PATH,
    FIGURES_PATH,
    tools,
    datasets,
    read_depths,
    mod1="bulk",
    mod2="single-cell",
):
    for tool in tools:
        for dataset in datasets:
            if dataset == "hgsoc-new":
                mod1 = "bulk_chunk_ribo"
                mod2 = "bulk_diss_polyA"
            elif dataset == "wilms_tumor":
                mod2 = "single-nucleus"
            for rd in read_depths:

                # Create plots
                bulk_vs_singlecell_matrix_viz(
                    DATA_PATH,
                    pseudobulk=False,
                    tool=tool,
                    dataset=dataset,
                    ncells="null",
                    rd=rd,
                    mod1=mod1,
                    mod2=mod2,
                    save_fig=True,
                    FIGURES_PATH=FIGURES_PATH,
                )


def super_plot_heatmaps_real_data(
    DATA_PATH,
    FIGURES_PATH,
    tool,
    dataset,
    read_depths,
    mod1="bulk",
    mod2="single-cell",
    save_fig=False,
):

    fig_name = f"superplot_{dataset}_{tool}_{mod1}_vs_{mod2}_similarity_matrix"
    super_extreme_point = 0
    all_matrices = []

    # Find all the data for the plots
    for rd in read_depths:
        matrix = load_heatmap_data(DATA_PATH, False, tool, dataset, "null", rd, mod1, mod2)
        expected_matches = load_expected_matches_real_data(DATA_PATH, dataset)
        if matrix is not None:
            ordered_matrix = order_matrix_by_expected_matches(
                matrix, expected_matches, tool, mod1, mod2
            )
            filtered_matrix = ordered_matrix.loc[
                [idx for idx in ordered_matrix.index if mod1 in idx],
                [col for col in ordered_matrix.columns if mod2 in col],
            ]
            # Store matrix for later use in super plot
            all_matrices.append({"rd": rd, "matrix": filtered_matrix})

            # Calculate the extreme point over all matrices to use the same color scale for all heatmaps in the super plot
            extreme_point = max(
                abs(filtered_matrix.min().min()), abs(filtered_matrix.max().max())
            )
            if extreme_point > super_extreme_point:
                super_extreme_point = extreme_point

    # Loop through the data again to create the super plot
    fig, axes = plt.subplots(
        1,
        len(read_depths),
        figsize=(3 * len(read_depths), 5),
        sharex=True,
        sharey=True,
    )
    fig.subplots_adjust(wspace=0.05, top=1)
    fig.suptitle(f"{tool} - {dataset} similarity matrices", fontsize=16)

    for rd in read_depths:
        print(f"Processing read depth {rd} for plotting...")
        if len(read_depths) == 1:
            ax = axes
        else:
            ax = axes[read_depths.index(rd)]
        matrix_list = [info["matrix"] for info in all_matrices if info["rd"] == rd]

        if not matrix_list or matrix_list[0] is None:
            # Skip this subplot if no data available
            ax.axis("off")
            ax.text(
                0.5,
                0.5,
                "No data",
                ha="center",
                va="center",
                transform=ax.transAxes,
                fontsize=12,
            )
            print(f"No data for read depth {rd}, skipping plot.")
            continue

        matrix = matrix_list[0]

        if tool == "CrosscheckFingerprints":
            cmap = plt.get_cmap("RdBu")
            cmap.set_bad(color="lightgrey")
            cax = ax.imshow(
                matrix,
                cmap=cmap,
                vmin=-super_extreme_point,
                vmax=super_extreme_point,
            )
        else:
            cmap = plt.get_cmap("Oranges")
            cmap.set_bad(color="lightgrey")
            cax = ax.imshow(matrix, cmap=cmap)
        ax.set_xticks(np.arange(len(matrix.columns)))
        ax.set_yticks(np.arange(len(matrix.index)))
        ax.set_xticklabels(matrix.columns, rotation=45, ha="right", fontsize=10)
        ax.set_yticklabels(matrix.index, fontsize=10)
        ax.set_title(f"read depth filter {rd}", pad=10, fontsize=12)

    # Shared colorbar for all subplots
    cbar = fig.colorbar(
        cax, ax=axes.ravel().tolist(), fraction=0.05, pad=0.04, shrink=0.7
    )

    # Save figure after all subplots are complete
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


def heatmap_plot_accuracy_metrics_pseudobulk(
    df, dataset, metric, save_fig=False, FIGURES_PATH=""
):

    fig_name = f"pseudobulk_{dataset}_{metric}_heatmap"

    fig, axes = plt.subplots(2, 2, figsize=(8, 10))
    fig.subplots_adjust(hspace=0.2)
    fig.suptitle(f"{metric} Heatmaps - {dataset} pseudobulks", fontsize=14, y=0.94)
    cmap = plt.get_cmap("Blues")
    cmap.set_bad(color="lightgrey")

    tools = ["CrosscheckFingerprints", "HYSYS", "NGSCheckmate", "Vireo"]
    for i, tool in enumerate(tools):
        sub_matrix = df[(df["dataset"] == dataset) & (df["tool"] == tool)].pivot_table(
            index=["read depth"], columns="ncells", values=metric
        )

        ax = axes[i // 2, i % 2]

        cax = ax.imshow(sub_matrix, cmap=cmap, vmin=0, vmax=1, aspect=1)

        # Add text annotations with F1 values
        for row in range(len(sub_matrix.index)):
            for col in range(len(sub_matrix.columns)):
                value = sub_matrix.iloc[row, col]
                if not np.isnan(value):
                    text_color = "white" if value > 0.5 else "black"
                    ax.text(
                        col,
                        row,
                        f"{value:.2f}",
                        ha="center",
                        va="center",
                        color=text_color,
                        fontsize=10,
                    )

        ax.set_xticks(np.arange(len(sub_matrix.columns)))
        ax.set_yticks(np.arange(len(sub_matrix.index)))

        # Format x-tick labels in scientific notation
        xticklabels = []
        for x in sub_matrix.columns:
            exponent = int(np.log10(x))
            mantissa = x / (10**exponent)
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


def heatmap_plot_accuracy_metrics_pseudobulk_rd0(
    df, dataset, save_fig=False, FIGURES_PATH=""
):
    metrics = ["fraction_inconclusive", "f1", "precision", "recall"]
    fig_name = f"pseudobulk_{dataset}_accuracy_metrics_heatmap_rd0"

    fig, ax = plt.subplots(2, 2, figsize=(8, 10))
    cmap = plt.get_cmap("Blues")
    cmap.set_bad(color="lightgrey")

    for i, metric in enumerate(metrics):
        matrix = df[(df["dataset"] == dataset) & (df["read depth"] == 0)].pivot_table(
            index=["tool"], columns="ncells", values=metric
        )

        ax[i // 2, i % 2].imshow(matrix, cmap=cmap, vmin=0, vmax=1, aspect=1)

        # Add text annotations with F1 values
        for row in range(len(matrix.index)):
            for col in range(len(matrix.columns)):
                value = matrix.iloc[row, col]
                if not np.isnan(value):
                    text_color = "white" if value > 0.5 else "black"
                    ax[i // 2, i % 2].text(
                        col,
                        row,
                        f"{value:.2f}",
                        ha="center",
                        va="center",
                        color=text_color,
                        fontsize=10,
                    )

        ax[i // 2, i % 2].set_title(f"{metric}")

        # Format x-tick labels in scientific notation
        xticklabels = []
        for x in matrix.columns:
            exponent = int(np.log10(x))
            mantissa = x / (10**exponent)
            if mantissa == 1.0:
                xticklabels.append(f"$10^{{{exponent}}}$")
            else:
                xticklabels.append(f"${int(mantissa)} \\times 10^{{{exponent}}}$")

        ax[i // 2, i % 2].set_yticks(np.arange(len(matrix.index)))
        ax[i // 2, i % 2].set_xticks(np.arange(len(matrix.columns)))

        # Only show x-label on bottom row
        if i // 2 == 1:
            ax[i // 2, i % 2].set_xlabel("# of cells")
            ax[i // 2, i % 2].set_xticklabels(xticklabels, ha="center")
        else:
            ax[i // 2, i % 2].tick_params(axis="x", labelbottom=False)

        # Only show y-label on left column
        if i % 2 == 0:
            ax[i // 2, i % 2].set_ylabel("Tool")
            ax[i // 2, i % 2].set_yticklabels(matrix.index)
        else:
            ax[i // 2, i % 2].tick_params(axis="y", labelleft=False)
    fig.suptitle(f"{dataset} pseudobulks")

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


def accuracy_metrics_errorbar_plot_missing_samples(
    DATA_PATH,
    experiment,
    dataset,
    ncells,
    rd,
    n_samples_removed,
    n_iterations,
    FIGURES_PATH="",
    save_fig=False,
):
    """
    Create errorbar plot showing accuracy metrics (fraction inconclusive, F1, precision, recall) for each tool,
    with different colors for different numbers of samples removed.
    Error bars represent standard deviation across iterations.
    Separate subplots for each metric.

    input:
        - DATA_PATH: path to data
        - experiment: type of experiment (missing_samples or double_samples)
        - dataset: name of dataset
        - ncells: number of cells in pseudobulk
        - rd: read depth filter cut-off
        - n_samples_removed: list of numbers of samples removed (e.g. [1, 5, 10])
        - n_iterations: number of iterations for each number of samples removed (e.g. 5)
        - FIGURES_PATH: path to save figures
        - save_fig: whether to save figure or just return it

    output:
        - errorbar plot is shown (or saved if save_fig=True)
    """
    fig_name = f"pseudobulk_{dataset}_accuracy_metrics_errorbar_plot_{experiment}_rd{rd}_ncells{ncells}"
    df = accuracy_metrics_averaged_over_iterations(
        DATA_PATH, n_iterations, dataset, ncells, rd, n_samples_removed, experiment
    )
    # Get unique tools and n_samples_removed values
    tools = df["tool"].unique()
    n_samples_vals = sorted(df["n_samples_removed"].unique())

    # Create color map for different n_samples_removed values
    colors = plt.cm.viridis(np.linspace(0, 1, len(n_samples_vals)))

    # Set up plot
    fig, axes = plt.subplots(2, 2, figsize=(10, 6), sharex=True, sharey=True)

    # Width for offsetting points within each tool group
    offset_width = 0.2
    offsets = np.linspace(-offset_width, offset_width, len(n_samples_vals))

    # Plot each n_samples_removed (or added) group
    for i, n_samp in enumerate(n_samples_vals):
        df_subset = df[df["n_samples_removed"] == n_samp]
        # Create x positions (tool indices + offset)
        x_pos = np.arange(len(tools)) + offsets[i]

        # Get y values and errors in the same order as tools
        for j, metric in enumerate(
            ["fraction_inconclusive", "f1", "precision", "recall"]
        ):
            if j == 0:
                ax = axes[0, 0]
                ax.set_title("Fraction Inconclusive", fontsize=12)
            elif j == 1:
                ax = axes[0, 1]
                ax.set_title("F1 Score", fontsize=12)
            elif j == 2:
                ax = axes[1, 0]
                ax.set_title("Precision", fontsize=12)
            elif j == 3:
                ax = axes[1, 1]
                ax.set_title("Recall", fontsize=12)
            y = [
                (
                    df_subset[df_subset["tool"] == t][f"av_{metric}"].values[0]
                    if len(df_subset[df_subset["tool"] == t]) > 0
                    else np.nan
                )
                for t in tools
            ]
            y_error = [
                (
                    df_subset[df_subset["tool"] == t][f"sd_{metric}"].values[0]
                    / np.sqrt(n_iterations)
                    if len(df_subset[df_subset["tool"] == t]) > 0
                    else np.nan
                )
                for t in tools
            ]
            # Plot with error bars
            ax.errorbar(
                x_pos,
                y,
                yerr=y_error,
                fmt="o",
                color=colors[i],
                ecolor=colors[i],
                elinewidth=2,
                capsize=4,
                alpha=0.7,
                label=f"{n_samp}",
            )

            # Customize plot
            ax.set_xticks(np.arange(len(tools)))
            ax.set_xticklabels(tools, rotation=0, ha="center")
            ax.set_ylim(0, 1.05)
            if j == 0:
                if experiment == "missing_samples":
                    ax.legend(title="# of samples removed", frameon=False)
                elif experiment == "double_samples":
                    ax.legend(title="# of double samples", frameon=False)

    plt.tight_layout()
    if save_fig:
        plt.savefig(
            os.path.join(
                FIGURES_PATH,
                f"{fig_name}.png",
            ),
            bbox_inches="tight",
            dpi=300,
        )
        plt.savefig(
            os.path.join(
                FIGURES_PATH,
                f"{fig_name}.pdf",
            ),
            bbox_inches="tight",
            dpi=300,
        )
    return


def heatmap_plot_accuracy_metrics_uneven_pseudobulk(
    df, dataset, metric, save_fig=False, FIGURES_PATH=""
):

    fig_name = f"pseudobulk_uneven_{dataset}_{metric}_heatmap"

    fig, axes = plt.subplots(2, 2, figsize=(10, 10))
    fig.subplots_adjust(hspace=0.2)
    fig.suptitle(f"{metric} heatmaps - {dataset} pseudobulks", fontsize=14, y=0.94)
    cmap = plt.get_cmap("Blues")
    cmap.set_bad(color="lightgrey")

    tools = ["CrosscheckFingerprints", "HYSYS", "NGSCheckmate", "Vireo"]
    for i, tool in enumerate(tools):
        sub_matrix = df[(df["dataset"] == dataset) & (df["tool"] == tool)].pivot_table(
            index=["read depth"], columns="ncells_2", values=metric
        )

        ax = axes[i // 2, i % 2]

        cax = ax.imshow(sub_matrix, cmap=cmap, vmin=0, vmax=1, aspect=1)

        # Add text annotations with F1 values
        for row in range(len(sub_matrix.index)):
            for col in range(len(sub_matrix.columns)):
                value = sub_matrix.iloc[row, col]
                if not np.isnan(value):
                    text_color = "white" if value > 0.5 else "black"
                    ax.text(
                        col,
                        row,
                        f"{value:.2f}",
                        ha="center",
                        va="center",
                        color=text_color,
                        fontsize=10,
                    )

        ax.set_xticks(np.arange(len(sub_matrix.columns)))
        ax.set_yticks(np.arange(len(sub_matrix.index)))

        # Format x-tick labels in scientific notation
        xticklabels = []
        for x in sub_matrix.columns:
            exponent = int(np.log10(x))
            mantissa = x / (10**exponent)
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


def super_plot_heatmaps_uneven_pseudobulk(
    DATA_PATH,
    FIGURES_PATH,
    tool,
    dataset,
    ncells_1,
    ncells_2_list,
    read_depths,
    save_fig=False,
):

    fig_name = f"superplot_pseudobulk_{dataset}_{tool}_all_samples_similarity_matrix"
    super_extreme_point = 0
    all_matrices = []

    # Find all the data for the plots
    for ncells_2 in ncells_2_list:
        for rd in read_depths:
            matrix = load_heatmap_data(
                DATA_PATH + f"uneven_pseudobulk_sizes/ncells_{ncells_2}/",
                True,
                tool,
                dataset,
                ncells_1,
                rd,
            )

            if matrix is not None:
                pseudobulk_1 = "_1"
                pseudobulk_2 = "_2"
                matrix_filtered = matrix.loc[
                    [idx for idx in matrix.index if idx.endswith(pseudobulk_1)],
                    [col for col in matrix.columns if col.endswith(pseudobulk_2)],
                ]

                # Store matrix for later use in super plot
                all_matrices.append(
                    {"ncells": ncells_2, "rd": rd, "matrix": matrix_filtered}
                )

                # Calculate the extreme point over all matrices to use the same color scale for all heatmaps in the super plot
                extreme_point = max(
                    abs(matrix_filtered.min().min()), abs(matrix_filtered.max().max())
                )
                if extreme_point > super_extreme_point:
                    super_extreme_point = extreme_point

    # Loop through the data again to create the super plot
    fig, axes = plt.subplots(
        len(read_depths),
        len(ncells_2_list),
        figsize=(5 * len(ncells_2_list), 4 * len(read_depths)),
        sharex=True,
        sharey=True,
    )
    fig.subplots_adjust(hspace=0.05, wspace=0.05, top=0.95)
    if len(read_depths) == 1:
        yval = 1.1
    else:
        yval = 0.98
    fig.suptitle(
        f"{tool} - {dataset} pseudobulks similarity matrices", fontsize=16, y=yval
    )

    for ncells in ncells_2_list:
        for rd in read_depths:
            if len(read_depths) == 1:
                ax = axes[ncells_2_list.index(ncells)]
            elif len(ncells_2_list) == 1:
                ax = axes[read_depths.index(rd)]
            else:
                ax = axes[read_depths.index(rd), ncells_2_list.index(ncells)]
            matrix_list = [
                info["matrix"]
                for info in all_matrices
                if info["ncells"] == ncells and info["rd"] == rd
            ]

            if not matrix_list or matrix_list[0] is None:
                # Skip this subplot if no data available
                ax.axis("off")
                ax.text(
                    0.5,
                    0.5,
                    "No data",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                    fontsize=12,
                )
                continue

            matrix = matrix_list[0]

            if tool == "CrosscheckFingerprints":
                cmap = plt.get_cmap("RdBu")
                cmap.set_bad(color="lightgrey")
                cax = ax.imshow(
                    matrix,
                    cmap=cmap,
                    vmin=-super_extreme_point,
                    vmax=super_extreme_point,
                )
            else:
                cmap = plt.get_cmap("Oranges")
                cmap.set_bad(color="lightgrey")
                cax = ax.imshow(matrix, cmap=cmap)
            ax.set_xticks(np.arange(len(matrix.columns)))
            ax.set_yticks(np.arange(len(matrix.index)))
            ax.set_xticklabels([""] * len(matrix.columns))
            ax.set_yticklabels([""] * len(matrix.index))
            # Add title with ncells to the top row
            if rd == read_depths[0]:
                ax.set_title(f"{ncells} cells", pad=10, fontsize=12)
            # Add text with read depth to the right column
            if ncells == ncells_2_list[-1]:
                ax.text(
                    1.05,
                    0.5,
                    f"Read depth filter: {rd}",
                    transform=ax.transAxes,
                    fontsize=12,
                    rotation=90,
                    va="center",
                )

    # Add one common colorbar for all subplots
    # cbar = fig.colorbar(
    #     cax, ax=axes.ravel().tolist(), fraction=0.046, pad=0.05, shrink=0.5
    # )

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


def heatmap_plot_accuracy_metrics_pseudobulk_vs_sc(
    df, dataset, metrics, save_fig=False, FIGURES_PATH=""
):

    fig_name = f"pseudobulk_vs_sc_{dataset}_accuracy_metrics_heatmap"

    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True, sharey=True)
    fig.subplots_adjust(hspace=0.2)
    fig.suptitle(
        f"Accuracy Metrics - {dataset} pseudobulk vs single-cell", fontsize=14, y=0.94
    )
    cmap = plt.get_cmap("Blues")
    cmap.set_bad(color="lightgrey")

    for i, metric in enumerate(metrics):
        sub_matrix = df[df["dataset"] == dataset].pivot_table(
            index=["tool"], columns="read depth", values=metric
        )

        ax = axes[i // 2, i % 2]

        cax = ax.imshow(sub_matrix, cmap=cmap, vmin=0, vmax=1, aspect=1)

        # Add text annotations with metric values
        for row in range(len(sub_matrix.index)):
            for col in range(len(sub_matrix.columns)):
                value = sub_matrix.iloc[row, col]
                if not np.isnan(value):
                    text_color = "white" if value > 0.5 else "black"
                    ax.text(
                        col,
                        row,
                        f"{value:.2f}",
                        ha="center",
                        va="center",
                        color=text_color,
                        fontsize=10,
                    )

        ax.set_xticks(np.arange(len(sub_matrix.columns)))
        ax.set_yticks(np.arange(len(sub_matrix.index)))
        ax.set_xticklabels(sub_matrix.columns, ha="center")
        ax.set_yticklabels(sub_matrix.index)

        # Only show x-label on bottom row
        if i // 2 == 1:
            ax.set_xlabel("Read depth cut-off")

        # Only show y-label on left column
        if i % 2 == 0:
            ax.set_ylabel("Tools")

        ax.set_title(f"{metric}", pad=10)

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
