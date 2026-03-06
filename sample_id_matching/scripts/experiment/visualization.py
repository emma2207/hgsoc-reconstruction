import os
import numpy as np
import matplotlib.pyplot as plt

from analysis_shared_functions import (
    parse_heatmap_matrix_bamixchecker,
    parse_heatmap_matrix_crosscheckfingerprints,
    parse_heatmap_matrix_hysys,
    parse_heatmap_matrix_ngscheckmate,
    parse_heatmap_matrix_vireo,
    load_expected_matches_real_data,
)


def load_heatmap_data(DATA_PATH, pseudobulk, tool, dataset, ncells, rd):

    print(f"Processing tool {tool}, dataset {dataset}, read depth {rd}...")
    try:
        # Create heatmap matrix
        if tool == "BAMixChecker":
            _, matrix = parse_heatmap_matrix_bamixchecker(
                DATA_PATH, pseudobulk, dataset, ncells
            )
        elif tool == "CrosscheckFingerprints":
            _, matrix = parse_heatmap_matrix_crosscheckfingerprints(
                DATA_PATH, pseudobulk, dataset, ncells, rd
            )
        elif tool == "HYSYS":
            _, matrix = parse_heatmap_matrix_hysys(
                DATA_PATH, pseudobulk, dataset, ncells, rd
            )
        elif tool == "NGSCheckmate":
            _, matrix = parse_heatmap_matrix_ngscheckmate(
                DATA_PATH, pseudobulk, dataset, ncells, rd
            )
        elif tool == "Vireo":
            matrix = parse_heatmap_matrix_vireo(
                DATA_PATH, pseudobulk, dataset, ncells, rd
            )
        else:
            print(f"Error! Do not recognize tool {tool}")
            print(
                "Choice of tool must be one of BAMixChecker, CrosscheckFingerprints, HYSYS, NGSCheckmate, or Vireo."
            )
            matrix = None
    except FileNotFoundError:
        print(
            f"Could not find data for {tool}, real dataset {dataset}, read depth {rd}. "
        )
        matrix = None

    return matrix


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
    matrix_df = load_heatmap_data(DATA_PATH, pseudobulk, tool, dataset, ncells, rd)
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
        plt.close()

    return


def all_samples_matrix_viz(
    DATA_PATH, pseudobulk, tool, dataset, ncells, rd, save_fig=False, FIGURES_PATH=""
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

    matrix_df = load_heatmap_data(DATA_PATH, pseudobulk, tool, dataset, ncells, rd)
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


def loop_heatmap_plots_pseudobulks(
    DATA_PATH, FIGURES_PATH, tools, datasets, ncells_list, read_depths
):

    for tool in tools:
        for dataset in datasets:
            for ncells in ncells_list:
                for rd in read_depths:
                    matrix = load_heatmap_data(
                        DATA_PATH, True, tool, dataset, ncells, rd
                    )

                    if matrix is not None:
                        pseudobulk_1 = "_1"
                        pseudobulk_2 = "_2"
                        matrix_filtered = matrix.loc[
                            [idx for idx in matrix.index if idx.endswith(pseudobulk_1)],
                            [
                                col
                                for col in matrix.columns
                                if col.endswith(pseudobulk_2)
                            ],
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
                            FIGURES_PATH=FIGURES_PATH,
                        )


def super_plot_heatmaps_pseudobulk(
    DATA_PATH, FIGURES_PATH, tool, dataset, ncells_list, read_depths, save_fig=False
):

    fig_name = f"superplot_pseudobulk_{dataset}_{tool}_all_samples_similarity_matrix"
    super_extreme_point = 0
    all_matrices = []

    # Find all the data for the plots
    for ncells in ncells_list:
        for rd in read_depths:
            matrix = load_heatmap_data(DATA_PATH, True, tool, dataset, ncells, rd)

            if matrix is not None:
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
    cbar = fig.colorbar(
        cax, ax=axes.ravel().tolist(), fraction=0.046, pad=0.05, shrink=0.5
    )

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


def order_matrix_by_expected_matches(
    matrix, expected_matches, tool, mod1="bulk", mod2="single-cell"
):  
    """
    Order the rows and columns of the matrix according to the expected matches between mod1 and mod2 samples, 
    so that samples that are expected to match show up on the diagonal of the heatmap.
    Samples with no expected match will be shown after the expected matches

    input:
        - matrix: dataframe containing similarity measures of all samples against all samples
        - expected_matches: dataframe containing expected matches between mod1 and mod2 samples, with one column for mod1 sample names and one column for mod2 sample names
        - tool: name of tool (used to determine how to interpret expected matches, e.g. Vireo has different expected match format than other tools)
        - mod1: name of first modality (e.g. "bulk")
        - mod2: name of second modality (e.g. "single-cell")

    output:
        - matrix ordered according to expected matches
    """
    # Parse expected matches to get lists of expected mod1 and mod2 samples
    mod1_col = [col for col in expected_matches.columns if mod1 in col]
    assert len(mod1_col) == 1, f"Expected exactly one {mod1} column in the matrix"
    mod2_col = [col for col in expected_matches.columns if mod2 in col]
    assert len(mod2_col) == 1, f"Expected exactly one {mod2} column in the matrix"

    mod1_samples = list(expected_matches[mod1_col[0]].dropna())
    mod2_samples = list(expected_matches[mod2_col[0]].dropna())
    matrix_samples = set(list(matrix.index) + list(matrix.columns))

    # Expected sample names should be prefixes of matrix sample names
    all_expected_samples = mod1_samples + mod2_samples
    all_expected_samples = [str(sample) for sample in all_expected_samples]
    print("All expected samples: " + str(all_expected_samples))

    # Find matching matrix samples for each expected sample
    # Iterate through expected samples in order to preserve ordering
    # Note: Not all expected samples may be present in the matrix
    ordered_matrix_samples = []
    for expected_sample in all_expected_samples:
        # Find matrix samples that contain this expected sample name as a substring
        matching_matrix_samples = [
            s for s in matrix_samples if str(expected_sample) in s
        ]
        if matching_matrix_samples:
            ordered_matrix_samples.extend(matching_matrix_samples)
        else:
            print(
                f"No matrix sample found for expected sample '{expected_sample}' (may not be in this batch)"
            )

    ordered_mod1_samples = [
        s
        for s in ordered_matrix_samples
        if any(str(mod1_sample) in s for mod1_sample in mod1_samples)
    ]
    ordered_mod2_samples = [
        s
        for s in ordered_matrix_samples
        if any(str(mod2_sample) in s for mod2_sample in mod2_samples)
    ]

    # Reorder matrix with matched samples
    if tool != "Vireo":
        matrix = matrix.loc[ordered_matrix_samples, ordered_matrix_samples]
    else:
        matrix = matrix.loc[ordered_mod1_samples, ordered_mod2_samples]

    return matrix


def loop_heatmap_plots_real_data(
    DATA_PATH,
    FIGURES_PATH,
    tools,
    datasets,
    read_depths,
    mod1="bulk",
    mod2="single-cell",
):
    ncells = "null"

    for tool in tools:
        for dataset in datasets:
            if dataset == "hgsoc-new":
                mod1 = "bulk_chunk_ribo"
                mod2 = "bulk_diss_polyA"
            elif dataset == "wilms_tumor":
                mod2 = "single-nucleus"
            for rd in read_depths:
                matrix = load_heatmap_data(DATA_PATH, False, tool, dataset, ncells, rd)

                # Load expected matches for real data
                expected_matches = load_expected_matches_real_data(DATA_PATH, dataset)
                ordered_matrix = order_matrix_by_expected_matches(
                    matrix, expected_matches, tool, mod1, mod2
                )

                # Create plots
                bulk_vs_singlecell_matrix_viz(
                    ordered_matrix,
                    pseudobulk=False,
                    tool=tool,
                    dataset=dataset,
                    ncells=ncells,
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
        matrix = load_heatmap_data(DATA_PATH, False, tool, dataset, "null", rd)
        expected_matches = load_expected_matches_real_data(DATA_PATH, dataset)
        if matrix is not None:
            ordered_matrix = order_matrix_by_expected_matches(
                matrix, expected_matches, tool, mod1, mod2
            )
            filtered_matrix = ordered_matrix.loc[
                [idx for idx in ordered_matrix.index if mod1 in idx],
                [col for col in ordered_matrix.columns if mod2 in col]
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

    fig, axes = plt.subplots(2, 2, figsize=(10, 10))
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

    fig, ax = plt.subplots(2, 2, figsize=(10, 10))
    cmap = plt.get_cmap("Blues")
    cmap.set_bad(color="lightgrey")

    for i, metric in enumerate(metrics):
        matrix = df[df["dataset"] == dataset].pivot_table(
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
