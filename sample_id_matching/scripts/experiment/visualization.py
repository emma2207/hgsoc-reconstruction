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
    sub_matrix = matrix_df[[col for col in matrix_df.columns if sc in col]]
    sub_matrix = sub_matrix.loc[
        [idx for idx in sub_matrix.index if "bulk" in idx and "diss" not in idx]
    ]

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
    ax.set_ylabel("bulk samples")
    ax.set_title(f"{tool} bulk vs {sc} similarity matrix", pad=20)

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
                    try:
                        # Create heatmap matrix
                        if tool == "BAMixChecker":
                            _, matrix = parse_heatmap_matrix_bamixchecker(
                                DATA_PATH, True, dataset, ncells
                            )
                        elif tool == "CrosscheckFingerprints":
                            _, matrix = parse_heatmap_matrix_crosscheckfingerprints(
                                DATA_PATH, True, dataset, ncells, rd
                            )
                        elif tool == "HYSYS":
                            _, matrix = parse_heatmap_matrix_hysys(
                                DATA_PATH, True, dataset, ncells, rd
                            )
                        elif tool == "NGSCheckmate":
                            _, matrix = parse_heatmap_matrix_ngscheckmate(
                                DATA_PATH, True, dataset, ncells, rd
                            )
                        elif tool == "Vireo":
                            matrix = parse_heatmap_matrix_vireo(
                                DATA_PATH, True, dataset, ncells, rd
                            )
                        else:
                            print(f"Error! Do not recognize tool {tool}.")
                            print(
                                "Choice of tool must be one of BAMixChecker, CrosscheckFingerprints, HYSYS, NGSCheckmate, or Vireo."
                            )
                            continue
                    except FileNotFoundError:
                        print(
                            f"Could not find data for {tool}, pseudobulk dataset {dataset}, ncells {ncells}, read depth {rd}. "
                        )
                        continue

                    pseudobulk_1 = "_1"
                    pseudobulk_2 = "_2"
                    matrix_filtered = matrix.loc[
                        [idx for idx in matrix.index if idx.endswith(pseudobulk_1)],
                        [col for col in matrix.columns if col.endswith(pseudobulk_2)],
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
            try:
                # Create heatmap matrix
                if tool == "BAMixChecker":
                    _, matrix = parse_heatmap_matrix_bamixchecker(
                        DATA_PATH, True, dataset, ncells
                    )
                elif tool == "CrosscheckFingerprints":
                    _, matrix = parse_heatmap_matrix_crosscheckfingerprints(
                        DATA_PATH, True, dataset, ncells, rd
                    )
                elif tool == "HYSYS":
                    _, matrix = parse_heatmap_matrix_hysys(
                        DATA_PATH, True, dataset, ncells, rd
                    )
                elif tool == "NGSCheckmate":
                    _, matrix = parse_heatmap_matrix_ngscheckmate(
                        DATA_PATH, True, dataset, ncells, rd
                    )
                elif tool == "Vireo":
                    matrix = parse_heatmap_matrix_vireo(
                        DATA_PATH, True, dataset, ncells, rd
                    )
                else:
                    print(f"Error! Do not recognize tool {tool}.")
                    print(
                        "Choice of tool must be one of BAMixChecker, CrosscheckFingerprints, HYSYS, NGSCheckmate, or Vireo."
                    )
                    continue
            except FileNotFoundError:
                print(
                    f"Could not find data for {tool}, pseudobulk dataset {dataset}, ncells {ncells}, read depth {rd}. "
                )
                matrix_filtered = None
                continue

            pseudobulk_1 = "_1"
            pseudobulk_2 = "_2"
            matrix_filtered = matrix.loc[
                [idx for idx in matrix.index if idx.endswith(pseudobulk_1)],
                [col for col in matrix.columns if col.endswith(pseudobulk_2)],
            ]

            # Store matrix for later use in super plot
            all_matrices.append({"ncells": ncells, "rd": rd, "matrix": matrix_filtered})

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
    fig.suptitle(
        f"{tool} - {dataset} pseudobulks similarity matrices", fontsize=16, y=0.98
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


def loop_heatmap_plots_real_data(
    DATA_PATH, FIGURES_PATH, tools, datasets, read_depths, mod1="bulk", mod2="single-cell"
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
                print(f"Processing tool {tool}, dataset {dataset}, read depth {rd}...")
                try:
                    # Create heatmap matrix
                    if tool == "BAMixChecker":
                        _, matrix = parse_heatmap_matrix_bamixchecker(
                            DATA_PATH, False, dataset, ncells
                        )
                    elif tool == "CrosscheckFingerprints":
                        _, matrix = parse_heatmap_matrix_crosscheckfingerprints(
                            DATA_PATH, False, dataset, ncells, rd
                        )
                    elif tool == "HYSYS":
                        _, matrix = parse_heatmap_matrix_hysys(
                            DATA_PATH, False, dataset, ncells, rd
                        )
                    elif tool == "NGSCheckmate":
                        _, matrix = parse_heatmap_matrix_ngscheckmate(
                            DATA_PATH, False, dataset, ncells, rd
                        )
                    elif tool == "Vireo":
                        matrix = parse_heatmap_matrix_vireo(
                            DATA_PATH, False, dataset, ncells, rd
                        )
                    else:
                        print(f"Error! Do not recognize tool {tool}")
                        print(
                            "Choice of tool must be one of BAMixChecker, CrosscheckFingerprints, HYSYS, NGSCheckmate, or Vireo."
                        )
                        continue
                except FileNotFoundError:
                    print(
                        f"Could not find data for {tool}, real dataset {dataset}, read depth {rd}. "
                    )
                    continue

                # Load expected matches for real data
                expected_matches = load_expected_matches_real_data(DATA_PATH, dataset)

                # Then use the ordering of the expected matches to order the matrix rows and columns before plotting,
                # so that samples that are expected to match show up on the diagonal
                # Samples with no expected match will be shown after the expected matches
                # bulk_col = [col for col in expected_matches.columns if bulk in col and "diss" not in col]
                # assert (
                #     len(bulk_col) == 1
                # ), "Expected exactly one bulk column in the matrix"
                # single_col = [
                #     col for col in expected_matches.columns if "single" in col or "diss" in col
                # ]
                # assert (
                #     len(single_col) == 1
                # ), "Expected exactly one single-cell/single-nucleus column in the matrix"

                # bulk_samples = list(expected_matches[bulk_col[0]].dropna())
                # single_samples = list(expected_matches[single_col[0]].dropna())

                # Filter matrix to only include samples matching expected samples
                # Expected sample names should be prefixes of matrix sample names
                # all_expected_samples = bulk_samples + single_samples
                # all_expected_samples = [str(sample) for sample in all_expected_samples]
                # matrix_samples = list(matrix.index)

                # Find matching matrix samples for each expected sample
                # ordered_matrix_samples = []
                # for expected_sample in all_expected_samples:
                #     if dataset != "hgsoc-new":
                #         matching_samples = [
                #             s for s in matrix_samples if expected_sample in s
                #         ]
                #     else:
                #         matching_bulk_samples = [
                #             s for s in bulk_samples if expected_sample in s
                #         ]
                #         matching_diss_bulk_samples = [
                #             s for s in single_samples if expected_sample in s
                #         ]
                #         matching_samples = matching_bulk_samples + matching_diss_bulk_samples
                #     if matching_samples:
                #         print(
                #             f"Found {len(matching_samples)} matches for {expected_sample}: {matching_samples}"
                #         )
                #         ordered_matrix_samples.extend(matching_samples)
                #     else:
                #         print(
                #             f"Warning: No match found in matrix for expected sample '{expected_sample}'"
                #         )

                # Reorder matrix with matched samples
                # matrix = matrix.loc[ordered_matrix_samples, ordered_matrix_samples]

                # Create plots
                bulk_vs_singlecell_matrix_viz(
                    matrix,
                    pseudobulk=False,
                    tool=tool,
                    dataset=dataset,
                    ncells=ncells,
                    rd=rd,
                    sc="bulk_diss",
                    save_fig=True,
                    FIGURES_PATH=FIGURES_PATH,
                )


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
