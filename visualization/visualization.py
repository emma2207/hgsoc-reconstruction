import os
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.cm import ScalarMappable

from analysis_shared_functions import (
    load_heatmap_data,
    order_matrix_by_expected_matches,
)
from accuracy_functions import (
    load_expected_matches_real_data,
    accuracy_metrics_averaged_over_iterations,
    count_matches_real_data,
    loop_accuracy_calculations,
)


METRIC_LABELS = {
    "fraction_inconclusive": "Fraction Inconclusive",
    "f1": "F1 Score",
    "precision": "Precision",
    "recall": "Recall",
}

TOOL_LABELS = {
    "Conpair": "Conpair",
    "CrosscheckFingerprints": "CrosscheckFingerprints",
    "HYSYS": "HYSYS",
    "NGSCheckmate": "NGSCheckMate",
    "ntsm": "ntsm",
    "OmicsPrint": "OmicsPrint",
    "Peddy": "Peddy",
    "Somalier": "Somalier",
    "TimeAttackGenComp": "TimeAttackGenComp",
    "Vireo": "Vireo",
}

DATASET_LABELS_SHORT = {
    "hgsoc": "HGSOC (Pilot)",
    "hgsoc-new": "HGSOC (New)",
    "high_grade_glioma": "HGG",
    "low_grade_glioma": "LGG",
    "wilms_tumor": "WT",
}

DATASET_LABELS = {
    "hgsoc": "High-Grade Serous Ovarian Cancer (Pilot)",
    "hgsoc-new": "High-Grade Serous Ovarian Cancer (New)",
    "high_grade_glioma": "High-Grade Glioma",
    "low_grade_glioma": "Low-Grade Glioma",
    "wilms_tumor": "Wilms Tumor",
}

MODALITY_LABELS = {
    "bulk": "Bulk",
    "single-cell": "Single-Cell",
    "single-nucleus": "Single-Nucleus",
    "bulk_chunk_ribo": "rRNA- Chunk Bulk",
    "bulk_diss_polyA": "Poly A+ Dissociated Bulk",
    "bulk_diss_ribo": "rRNA- Dissociated Bulk",
    "bulk_dissociated_polyA": "Poly A+ Dissociated Bulk",
    "bulk_dissociated_ribo": "rRNA- Dissociated Bulk",
}


def _build_real_data_rd_tag(mod1, mod2, rd, rd_mod1, rd_mod2):
    """Build read-depth selector passed to data loaders for real-data runs."""
    if rd_mod1 is not None and rd_mod2 is not None:
        return (
            f"{mod1}_{rd_mod1}_{mod2}_{rd_mod2}",
            f"{mod1}:{rd_mod1}, {mod2}:{rd_mod2}",
        )
    return rd, str(rd)


def _modality_in_label_local(modality, label):
    modality_key = str(modality).lower().replace("-", "_").replace(" ", "_")
    label_key = str(label).lower().replace("-", "_").replace(" ", "_")

    if modality_key in label_key:
        return True

    modality_tokens = modality_key.split("_")
    if "single" in modality_tokens and "single" in label_key:
        return True

    return False


def _modality_submatrix(matrix, mod1, mod2):
    """Return matrix restricted to mod1 rows and mod2 columns, correcting swapped axes."""

    def normalize_key(value):
        return re.sub(r"[^a-z0-9]+", "_", str(value).lower()).strip("_")

    def has_modality(label, modality_key):
        label_key = normalize_key(label)
        if modality_key in label_key:
            return True
        # Pooled single-cell labels are sometimes represented as "poolX" only.
        if "pooled" in modality_key and label_key.startswith("pool"):
            return True
        return False

    mod1_key = normalize_key(mod1)
    mod2_key = normalize_key(mod2)

    row_has_mod1 = sum(has_modality(label, mod1_key) for label in matrix.index)
    row_has_mod2 = sum(has_modality(label, mod2_key) for label in matrix.index)
    col_has_mod1 = sum(has_modality(label, mod1_key) for label in matrix.columns)
    col_has_mod2 = sum(has_modality(label, mod2_key) for label in matrix.columns)

    # Some outputs can be transposed depending on tool/parser behavior.
    if (
        row_has_mod1 == 0
        and col_has_mod1 > 0
        and row_has_mod2 > 0
        and col_has_mod2 == 0
    ):
        matrix = matrix.T

    mod1_rows = [label for label in matrix.index if has_modality(label, mod1_key)]
    mod2_cols = [label for label in matrix.columns if has_modality(label, mod2_key)]

    # If one modality is not explicitly encoded in labels, use the complement of
    # the other modality as a fallback partition for mixed-modality matrices.
    if len(mod2_cols) == 0 and len(mod1_rows) > 0:
        mod1_cols = [label for label in matrix.columns if has_modality(label, mod1_key)]
        mod2_cols = [label for label in matrix.columns if label not in mod1_cols]
    if len(mod1_rows) == 0 and len(mod2_cols) > 0:
        mod2_rows = [label for label in matrix.index if has_modality(label, mod2_key)]
        mod1_rows = [label for label in matrix.index if label not in mod2_rows]

    if len(mod1_rows) == 0 or len(mod2_cols) == 0:
        return matrix

    return matrix.loc[mod1_rows, mod2_cols]


def _subset_matrix_by_modalities(matrix, mod1, mod2):
    row_labels = [label for label in matrix.index if _modality_in_label_local(mod1, label)]
    col_labels = [
        label for label in matrix.columns if _modality_in_label_local(mod2, label)
    ]

    if len(row_labels) == 0 or len(col_labels) == 0:
        return matrix

    return matrix.loc[row_labels, col_labels]


def _tool_colormap(tool):
    cmap_name = "RdBu" if tool == "CrosscheckFingerprints" else "Oranges"
    cmap = plt.get_cmap(cmap_name).copy()
    cmap.set_bad(color="lightgrey")
    return cmap


def _tool_norm(tool, limit):
    if tool == "CrosscheckFingerprints":
        return colors.Normalize(vmin=-limit, vmax=limit)
    return colors.Normalize(vmin=0, vmax=limit)


def _plot_tool_heatmap(ax, matrix, tool, method="imshow", symmetric_limit=None):
    cmap = _tool_colormap(tool)
    plotter = getattr(ax, method)
    if tool == "CrosscheckFingerprints" and symmetric_limit is not None:
        return plotter(matrix, cmap=cmap, vmin=-symmetric_limit, vmax=symmetric_limit)
    return plotter(matrix, cmap=cmap)


def bulk_vs_singlecell_matrix_viz(
    DATA_PATH,
    pseudobulk,
    tool,
    dataset,
    ncells,
    rd=None,
    rd_mod1=None,
    rd_mod2=None,
    mod1="bulk",
    mod2="single-cell",
    save_fig=False,
    remove_missing_data=False,
    FIGURES_PATH="",
):
    """
    Docstring for bulk_vs_singlecell_matrix_viz

    :param matrix_df: Dataframe containing a heatmap / matrix of similarity measures of all samples against all samples
    :param pseudobulk: True / False
    :param tool: BAMixChecker / CrosscheckFingerprints / HYSYS / NGSCheckmate / Vireo
    :param dataset: dataset name
    :param ncells: number of cells in pseudobulk, "null" in real data
    :param rd: read depth filter cut-off (legacy / shared for both modalities)
    :param rd_mod1: read depth filter for first modality in real-data runs
    :param rd_mod2: read depth filter for second modality in real-data runs
    :param mod1: first modality (e.g., "bulk")
    :param mod2: second modality (e.g., "single-cell")
    :param save_fig: True / False
    :param FIGURES_PATH: path that figures get saved to
    :param remove_missing_data: True / False
    """

    fontsize=12
    rd_to_load, rd_label = _build_real_data_rd_tag(mod1, mod2, rd, rd_mod1, rd_mod2)

    if pseudobulk:
        fig_name = f"pseudobulk_{dataset}_{tool}_rd{rd}_ncells{ncells}_{mod1}_vs_{mod2}_similarity_matrix"
    else:
        fig_name = f"{dataset}_{tool}_rd{rd_to_load}_similarity_matrix"

    # Visualize bulk against single-cell or single-nucleus samples
    matrix_df = load_heatmap_data(
        DATA_PATH, pseudobulk, tool, dataset, ncells, rd_to_load, mod1, mod2
    )
    if matrix_df is None:
        print(
            f"No matrix data available for {tool} on {dataset} with read depth {rd_label} and ncells {ncells}. Skipping plot."
        )
        return

    sub_matrix = _modality_submatrix(matrix_df, mod1, mod2)
    if not pseudobulk:
        expected_matches = load_expected_matches_real_data(DATA_PATH, dataset)
        try:
            sub_matrix = order_matrix_by_expected_matches(
                sub_matrix,
                expected_matches,
                mod1,
                mod2,
            )
        except ValueError as exc:
            print(
                "Could not order by expected matches for requested modalities "
                f"({mod1} vs {mod2}): {exc}. Falling back to direct modality filtering."
            )
            sub_matrix = _subset_matrix_by_modalities(sub_matrix, mod1, mod2)

    extreme_point = max(abs(sub_matrix.min().min()), abs(sub_matrix.max().max()))

    if remove_missing_data:
        sub_matrix = sub_matrix.dropna(axis=0, how="all").dropna(axis=1, how="all")

    fig, ax = plt.subplots(figsize=(6, 6))
    cax = _plot_tool_heatmap(
        ax,
        sub_matrix,
        tool,
        method="matshow",
        symmetric_limit=extreme_point,
    )

    ax.set_xticks(np.arange(len(sub_matrix.columns)))
    ax.set_yticks(np.arange(len(sub_matrix.index)))
    # ax.set_xticklabels([""] * len(sub_matrix.columns))
    # ax.set_yticklabels([""] * len(sub_matrix.index))
    xlabels = [
        f"{col.replace('HGSOC-', '').replace(f'_{mod2}', '')}"
        for col in sub_matrix.columns
    ]
    ylabels = [
        f"{idx.replace('HGSOC-', '').replace(f'_{mod1}', '')}"
        for idx in sub_matrix.index
    ]
    ax.set_xticklabels(xlabels, rotation=90, ha="center")
    ax.set_yticklabels(ylabels)
    ax.xaxis.set_ticks_position("bottom")
    ax.xaxis.set_label_position("bottom")
    ax.set_xlabel(f"{MODALITY_LABELS.get(mod2, mod2)}", fontsize=fontsize)
    ax.set_ylabel(f"{MODALITY_LABELS.get(mod1, mod1)}", fontsize=fontsize)
    ax.set_title(f"HGSOC - {tool} Similarity Matrix", pad=20, fontsize=fontsize+2)

    # fig.colorbar(cax, fraction=0.046, pad=0.04, shrink=0.6)

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

    plt.close(fig)
    return


def all_samples_matrix_viz(
    DATA_PATH,
    pseudobulk,
    tool,
    dataset,
    ncells,
    rd,
    save_fig=False,
    FIGURES_PATH="",
    mod1="bulk",
    mod2="single-cell",
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

    matrix_df = load_heatmap_data(
        DATA_PATH, pseudobulk, tool, dataset, ncells, rd, mod1, mod2
    )
    if matrix_df is not None:
        # Visualize ALL samples against ALL samples
        extreme_point = max(abs(matrix_df.min().min()), abs(matrix_df.max().max()))

        fig, ax = plt.subplots(figsize=(10, 10))
        cax = _plot_tool_heatmap(
            ax,
            matrix_df,
            tool,
            method="matshow",
            symmetric_limit=extreme_point,
        )

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
    tools,
    datasets,
    ncells,
    rd,
    experiment,
    save_fig=False,
):

    fig_name = f"superplot_pseudobulk_ncells_{ncells}_rd_{rd}_similarity_matrix"
    fontsize = 12
    all_matrices = []

    # Find all the data for the plots
    for dataset in datasets:
        for tool in tools:
            crosscheck_extreme_point = None
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
                if tool == "CrosscheckFingerprints":
                    crosscheck_extreme_point = max(
                        abs(matrix_filtered.min().min()), abs(matrix_filtered.max().max())
                    )
                # Store matrix for later use in super plot
                all_matrices.append(
                    {
                        "tool": tool, 
                        "dataset": dataset, 
                        "matrix": matrix_filtered, 
                        "crosscheck_extreme_point": crosscheck_extreme_point,
                     }
                )


    # Loop through the data again to create the super plot
    fig, axes = plt.subplots(
        len(tools),
        len(datasets),
        figsize=(3 * len(datasets), 3 * len(tools)),
    )
    fig.subplots_adjust(hspace=0.05, wspace=0.05, top=0.95)
    if len(tools) == 1:
        yval = 1.1
    else:
        yval = 1.
    fig.suptitle(
        f"Pseudobulk Sample Similarity Matrices", fontsize=fontsize+2, y=yval
    )

    for dataset in datasets:
        for tool in tools:
            if len(tools) == 1:
                ax = axes[datasets.index(dataset)]
            elif len(datasets) == 1:
                ax = axes[tools.index(tool)]
            else:
                ax = axes[tools.index(tool), datasets.index(dataset)]
            matrix_list = [
                info["matrix"]
                for info in all_matrices
                if info["tool"] == tool and info["dataset"] == dataset
            ]
            crosscheck_extreme_point = [info["crosscheck_extreme_point"] for info in all_matrices if info["tool"] == tool and info["dataset"] == dataset][0]

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

            cax = _plot_tool_heatmap(
                ax,
                matrix,
                tool,
                method="imshow",
                symmetric_limit=crosscheck_extreme_point if tool == "CrosscheckFingerprints" else None
            )
            ax.set_xticks(np.arange(len(matrix.columns)))
            ax.set_yticks(np.arange(len(matrix.index)))
            ax.set_xticklabels([""] * len(matrix.columns))
            ax.set_yticklabels([""] * len(matrix.index))
            # Add title with ncells to the top row
            if tool == tools[0]:
                ax.set_title(f"{DATASET_LABELS[dataset]}", pad=10, fontsize=fontsize)
            # Add text with read depth to the left column
            if dataset == datasets[0]:
                ax.set_ylabel(f"{TOOL_LABELS[tool]}", fontsize=fontsize, rotation=90, labelpad=10)

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


def super_plot_heatmaps_fixed_rd_pseudobulk(
    DATA_PATH,
    FIGURES_PATH,
    tools,
    dataset,
    ncells_list,
    rd,
    experiment,
    save_fig=False,
):

    fig_name = f"superplot_pseudobulk_{dataset}_rd_{rd}_similarity_matrix_accuracy"
    fontsize = 20
    super_extreme_point = -np.inf
    tool_extreme_points = {t: 0 for t in tools}
    row_has_data = {t: False for t in tools}
    row_mappables: dict[str, ScalarMappable | None] = {t: None for t in tools}
    all_matrices = []

    # Find all the sample similarity data for the plots
    for ncells in ncells_list:
        for tool in tools:
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
                if matrix_filtered is not None and not matrix_filtered.empty:
                    all_matrices.append(
                        {"ncells": ncells, "tool": tool, "matrix": matrix_filtered}
                    )

                # in the "Find all the data for the plots" loop, keep current logic and also update per-tool max
                extreme_point = max(
                    abs(matrix_filtered.min().min()), abs(matrix_filtered.max().max())
                )
                if extreme_point > super_extreme_point:
                    super_extreme_point = extreme_point
                if extreme_point > tool_extreme_points[tool]:
                    tool_extreme_points[tool] = extreme_point

    # Calculate accuracy metrics for all heatmaps in the grid
    accuracy_df = loop_accuracy_calculations(
        DATA_PATH=DATA_PATH,
        experiment="pseudobulk",
        tools=tools,
        datasets=[dataset],
        ncells_list=ncells_list,
        read_depths=[rd],
    )

    # Loop through the data again to create the super plot
    fig, axes = plt.subplots(
        len(tools),
        len(ncells_list),
        figsize=(5 * len(ncells_list), 4 * len(tools)),
        sharex=True,
        sharey=True,
    )
    fig.subplots_adjust(hspace=0.05, wspace=0.05, top=0.95)
    if len(tools) == 1:
        yval = 1.1
    else:
        yval = 1.0
    fig.suptitle(
        f"{dataset.replace('_', ' ')} pseudobulks - sample similarity".title(),
        fontsize=fontsize + 4,
        y=yval,
    )

    for ncells in ncells_list:
        for tool in tools:
            # Deal with the axes
            tool_mappable = None
            if len(tools) == 1 and len(ncells_list) == 1:
                ax = axes
                row_left_ax = ax
            elif len(tools) == 1:
                ax = axes[ncells_list.index(ncells)]
                row_left_ax = axes[0]
            elif len(ncells_list) == 1:
                ax = axes[tools.index(tool)]
                row_left_ax = axes[tools.index(tool)]
            else:
                ax = axes[tools.index(tool), ncells_list.index(ncells)]
                row_left_ax = axes[tools.index(tool), 0]
            matrix_list = [
                info["matrix"]
                for info in all_matrices
                if info["ncells"] == ncells and info["tool"] == tool
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
                    fontsize=fontsize,
                )
                continue

            matrix = matrix_list[0]

            # in plotting loop, replace imshow branch with row-normalized plotting
            cmap = _tool_colormap(tool)
            norm = _tool_norm(tool, tool_extreme_points[tool])

            ax.imshow(matrix, cmap=cmap, norm=norm)
            row_has_data[tool] = True
            tool_mappable = ScalarMappable(norm=norm, cmap=cmap)
            tool_mappable.set_array([])
            row_mappables[tool] = tool_mappable
            ax.set_xticks(np.arange(len(matrix.columns)))
            ax.set_yticks(np.arange(len(matrix.index)))
            ax.set_xticklabels([""] * len(matrix.columns))
            ax.set_yticklabels([""] * len(matrix.index))
            if ncells == ncells_list[0]:
                ax.set_ylabel(f"{TOOL_LABELS[tool]}", fontsize=fontsize, rotation=90, labelpad=10)
            # Add title with ncells to the top row
            if tool == tools[0]:
                ax.set_title(f"{ncells:,.0f} UMIs", pad=10, fontsize=fontsize)


            # Add text with accuracy metrics to each subplot
            accuracy_metrics = accuracy_df[(accuracy_df["tool"] == tool) & (accuracy_df["ncells"] == ncells)]
            if not accuracy_metrics.empty:
                ax.text(
                    0.03,
                    0.03,
                    f"F1: {accuracy_metrics['f1'].values[0]:.2f}\nRecall: {accuracy_metrics['recall'].values[0]:.2f}\nPrecision: {accuracy_metrics['precision'].values[0]:.2f}\n% Inconclusive: {accuracy_metrics['fraction_inconclusive'].values[0]:.2f}",
                    transform=ax.transAxes,
                    fontsize=fontsize - 2,
                    verticalalignment="bottom",
                    horizontalalignment="left",
                )

    # Add text with tool name to the left of each row
    # Needed if there are missing data for some tools, so that the row labels are still aligned with the correct rows
    # for itool, tool in enumerate(reversed(tools)):
    #     row_bbox = row_left_ax.get_position()
    #     fig.text(
    #         row_bbox.x0 - 0.02,
    #         row_bbox.y0 * (itool * 1.93 + 1) + row_bbox.height / 2,
    #         f"{tool}",
    #         transform=fig.transFigure,
    #         fontsize=fontsize,
    #         rotation=90,
    #         va="center",
    #         ha="right",
    #     )

    # one colorbar per row/tool
    for i, tool in enumerate(tools):
        if not row_has_data[tool] or row_mappables[tool] is None:
            continue

        if len(tools) == 1 and len(ncells_list) == 1:
            row_axes = [axes]
        elif len(tools) == 1:
            row_axes = axes  # 1D over columns
        elif len(ncells_list) == 1:
            row_axes = [axes[i]]  # 1D over rows
        else:
            row_axes = axes[i, :]  # full row

        row_mappable = row_mappables[tool]
        if row_mappable is None:
            continue

        cbar = fig.colorbar(
            row_mappable,
            ax=row_axes,
            fraction=0.08,
            pad=0.03,
            shrink=0.75,
            aspect=18,
        )
        cbar.ax.tick_params(labelsize=fontsize)

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

    plt.close(fig)

    return


def loop_heatmap_plots_real_data(
    DATA_PATH,
    FIGURES_PATH,
    tools,
    datasets,
    read_depths_mod1,
    read_depths_mod2,
    mod1="bulk",
    mod2="single-cell",
):
    if len(read_depths_mod1) != len(read_depths_mod2):
        raise ValueError(
            "read_depths_mod1 and read_depths_mod2 must have the same length."
        )

    original_mod1 = mod1
    original_mod2 = mod2
    for tool in tools:
        for dataset in datasets:
            if dataset == "hgsoc-new":
                mod1 = "bulk_chunk_ribo"
                mod2 = "bulk_diss_polyA"
            elif dataset == "wilms_tumor":
                mod2 = "single-nucleus"
            else:
                mod1 = original_mod1
                mod2 = original_mod2
            rd_pairs = list(zip(read_depths_mod1, read_depths_mod2))

            for rd1, rd2 in rd_pairs:

                # Create plots
                bulk_vs_singlecell_matrix_viz(
                    DATA_PATH,
                    pseudobulk=False,
                    tool=tool,
                    dataset=dataset,
                    ncells="null",
                    rd_mod1=rd1,
                    rd_mod2=rd2,
                    mod1=mod1,
                    mod2=mod2,
                    save_fig=True,
                    FIGURES_PATH=FIGURES_PATH,
                )


def super_plot_heatmaps_real_data(
    DATA_PATH,
    FIGURES_PATH,
    tools,
    dataset,
    read_depths_mod1,
    read_depths_mod2,
    mod1="bulk",
    mod2="single-cell",
    remove_missing_data=False,
    save_fig=False,
):
    if len(read_depths_mod1) != len(read_depths_mod2):
        raise ValueError(
            "read_depths_mod1 and read_depths_mod2 must have the same length."
        )

    fig_name = f"superplot_{dataset}_{mod1}_vs_{mod2}_similarity_matrix"
    fontsize = 12
    tool_extreme_points = np.full(len(tools), -np.inf)
    all_matrices = []
    row_mappables: dict[str, ScalarMappable | None] = {t: None for t in tools}
    rd_pairs = [
        (rd1, rd2, f"{mod1}_{rd1}_{mod2}_{rd2}")
        for rd1, rd2 in zip(read_depths_mod1, read_depths_mod2)
    ]

    # Find all the data for the plots
    for i, tool in enumerate(tools):
        for rd1, rd2, rd_to_load in rd_pairs:
            matrix = load_heatmap_data(
                DATA_PATH, False, tool, dataset, "null", rd_to_load, mod1, mod2
            )
            if matrix is not None:
                matrix = _modality_submatrix(matrix, mod1, mod2)
                expected_matches = load_expected_matches_real_data(DATA_PATH, dataset)
                try:
                    ordered_matrix = order_matrix_by_expected_matches(
                        matrix, expected_matches, mod1, mod2
                    )
                except ValueError as exc:
                    print(
                        "Could not order by expected matches for requested modalities "
                        f"({mod1} vs {mod2}): {exc}. Falling back to direct modality filtering."
                    )
                    ordered_matrix = _subset_matrix_by_modalities(matrix, mod1, mod2)
                if remove_missing_data:
                    ordered_matrix = ordered_matrix.dropna(axis=0, how="all")

                all_matrices.append(
                    {
                        "rd_mod1": rd1,
                        "rd_mod2": rd2,
                        "rd": rd_to_load,
                        "tool": tool,
                        "matrix": ordered_matrix,
                    }
                )

                extreme_point = max(
                    abs(ordered_matrix.min().min()), abs(ordered_matrix.max().max())
                )
                if extreme_point > tool_extreme_points[i]:
                    tool_extreme_points[i] = extreme_point

    fig, axes = plt.subplots(
        len(tools),
        len(rd_pairs),
        figsize=(2.5 * len(rd_pairs), 2 * len(tools)),
        sharex="col",
        sharey="row",
    )
    fig.subplots_adjust(wspace=0.10, hspace=0.1)
    fig.suptitle(f"{dataset.replace('_', ' ')} similarity matrices".title(), fontsize=fontsize + 4, y=0.95)
    # fig.suptitle("Wilms Tumor Similarity Matrices", fontsize=16, y=0.93)

    for i, tool in enumerate(tools):
        row_has_data = False
        tool_mappable = None

        cmap = _tool_colormap(tool)
        norm = _tool_norm(tool, tool_extreme_points[i])
        

        for rd1, rd2, rd_to_load in rd_pairs:
            if rd1 is not None and rd2 is not None:
                rd_label = f"{rd1}, {rd2}"
            else:
                rd_label = f"{rd_to_load}"
            print(f"Processing read depth {rd_label} and tool {tool} for plotting...")
            if len(rd_pairs) == 1 and len(tools) == 1:
                ax = axes
            elif len(rd_pairs) == 1:
                ax = axes[tools.index(tool)]
            elif len(tools) == 1:
                ax = axes[rd_pairs.index((rd1, rd2, rd_to_load))]
            else:
                ax = axes[tools.index(tool), rd_pairs.index((rd1, rd2, rd_to_load))]

            matrix_list = [
                info["matrix"]
                for info in all_matrices
                if info["rd"] == rd_to_load and info["tool"] == tool
            ]

            if not matrix_list or matrix_list[0] is None:
                ax.axis("off")
                ax.text(
                    0.5,
                    0.5,
                    "No data",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                    fontsize=fontsize,
                )
                print(
                    f"No data for read depth {rd_label} and tool {tool}, skipping plot."
                )
                continue

            matrix = matrix_list[0]
            row_has_data = True

            ax.imshow(matrix, cmap=cmap, norm=norm)
            tool_mappable = ScalarMappable(norm=norm, cmap=cmap)
            tool_mappable.set_array([])
            row_mappables[tool] = tool_mappable

            ax.set_xticks(np.arange(len(matrix.columns)))
            ax.set_yticks(np.arange(len(matrix.index)))
            if dataset == "hgsoc":
                xlabels = [x.replace(mod2, "") + f"{mod2}" for x in matrix.columns]
                ylabels = [y.replace(mod1, "") + f"{mod1}" for y in matrix.index]
            else:
                xlabels = matrix.columns
                xlabels = [x.replace(mod2, "").rstrip("_").replace("_missing", "") for x in xlabels]
                ylabels = matrix.index
                ylabels = [y.replace(mod1, "").rstrip("_") for y in ylabels]
            ax.set_xticklabels(xlabels, rotation=90, ha="right", fontsize=4)
            ax.set_yticklabels(ylabels, fontsize=4)
            ax.set_xticklabels([""] * len(matrix.columns))
            ax.set_yticklabels([""] * len(matrix.index))
            if tool == tools[0]:
                ax.yaxis.set_label_position("left")
                ax.set_title(
                    f"Read depth {rd_label}", fontsize=fontsize, pad=10
                )
            if (rd1, rd2, rd_to_load) == rd_pairs[0]:
                ax.set_ylabel(f"{tool}", fontsize=fontsize, rotation=90, labelpad=10, va="center")

        if row_has_data:
            if len(rd_pairs) == 1 and len(tools) == 1:
                row_axes = [axes]
            elif len(tools) == 1:
                row_axes = axes
            elif len(rd_pairs) == 1:
                row_axes = [axes[i]]
            else:
                row_axes = axes[i, :]

            if tool_mappable is not None:
                cbar = fig.colorbar(
                    tool_mappable,
                    ax=row_axes,
                    orientation="vertical",
                    pad=0.03,
                    fraction=0.08,
                    shrink=0.75,
                )
                cbar.ax.tick_params(labelsize=fontsize)

    if save_fig:
        fig.savefig(
            os.path.join(FIGURES_PATH, f"{fig_name}.png"),
            bbox_inches="tight",
            dpi=300,
        )
        fig.savefig(
            os.path.join(FIGURES_PATH, f"{fig_name}.pdf"),
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


def super_plot_heatmaps_real_data_multimodal_hgsoc(
    DATA_PATH,
    FIGURES_PATH,
    tool,
    rd,
    save_fig=False,
):

    fig_name = f"superplot_hgsoc_multimodal_similarity_matrix_{tool}_rd{rd}"
    mod_list1 = ["single-cell", "bulk_dissociated_polyA", "bulk_dissociated_ribo"] # columns
    mod_list2 = ["bulk_dissociated_polyA", "bulk_dissociated_ribo", "bulk_chunk_ribo"] # rows
    dataset = "hgsoc"
    all_matrices = []
    overall_extreme_point = 0
    col_mappable = None

    # Find all the data for the plots
    for i, mod1 in enumerate(mod_list1):
        for j, mod2 in enumerate(mod_list2):

            matrix = load_heatmap_data(
                DATA_PATH, False, tool, dataset, "null", rd, mod1, mod2
            )
            if matrix is None:
                # Try loading the transposed matrix in case of swapped axes in the output
                matrix = load_heatmap_data(
                    DATA_PATH, False, tool, dataset, "null", rd, mod2, mod1
                )

            if matrix is not None:
                print(f"Processing modalities {mod1} vs {mod2} for plotting...")
                # Clean up labels
                matrix.index = [
                    idx.replace(f"{mod1}_vs_{mod2}/", "")
                    .replace(f"{mod2}_vs_{mod1}/", "")
                    .replace(f"_{mod1}_{mod1}", f"_{mod1}")
                    .replace(f"_{mod2}_{mod2}", f"_{mod2}")
                    for idx in matrix.index
                ]
                matrix.columns = [
                    col.replace(f"{mod1}_vs_{mod2}/", "")
                    .replace(f"{mod2}_vs_{mod1}/", "")
                    .replace(f"_{mod1}_{mod1}", f"_{mod1}")
                    .replace(f"_{mod2}_{mod2}", f"_{mod2}")
                    for col in matrix.columns
                ]
                # Keep grid semantics consistent: rows correspond to mod_list2,
                # columns correspond to mod_list1.
                matrix = _modality_submatrix(matrix, mod2, mod1)

                all_matrices.append(
                    {
                        "mod1": mod1,
                        "mod2": mod2,
                        "matrix": matrix,
                    }
                )
                extreme_point = max(abs(matrix.min().min()), abs(matrix.max().max()))
                if extreme_point > overall_extreme_point:
                    overall_extreme_point = extreme_point

    # Plot data
    fig, axes = plt.subplots(3, 3, figsize=(8, 8), sharex=True, sharey=True)

    cmap = _tool_colormap(tool)
    norm = _tool_norm(tool, overall_extreme_point)

    for j, mod1 in enumerate(mod_list1):
        for i, mod2 in enumerate(mod_list2):
            if i < j:
                axes[i, j].axis("off")
                continue
            ax = axes[i, j]
            matrix_list = [
                info["matrix"]
                for info in all_matrices
                if info["mod1"] == mod1 and info["mod2"] == mod2
            ]

            if not matrix_list or matrix_list[0] is None:
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

            ax.imshow(matrix, cmap=cmap, norm=norm)
            col_mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
            col_mappable.set_array([])

            # Y-label on left column
            if j == 0:
                ax.set_ylabel(f"{MODALITY_LABELS.get(mod_list2[i], mod_list2[i])}", fontsize=12)

            # Modality on bottom row
            if i == len(mod_list2) - 1:
                ax.set_xlabel(f"{MODALITY_LABELS.get(mod_list1[j], mod_list1[j])}", fontsize=12)
            
            # Only show x-tick labels on bottom row
            if i == len(mod_list2) - 1:
                ax.set_xticks(np.arange(len(matrix.columns)))
                xlabels = [col.replace(f"_{mod1}", "").replace(f"_{mod2}", "") for col in matrix.columns]
                ax.set_xticklabels(xlabels, rotation=90)
            else:
                ax.set_xticks([])
            # Only show y-tick labels on left column
            if j == 2:
                ax.set_yticks(np.arange(len(matrix.index)))
                ylabels = [idx.replace(f"_{mod1}", "").replace(f"_{mod2}", "") for idx in matrix.index]
                ax.set_yticklabels(ylabels)
            else:
                ax.set_yticks([])
            ax.tick_params(axis='both', which='major', labelsize=10)
    
    # Add one common colorbar for all subplots
    cax = ax.inset_axes([0.3, 1.6, 0.15, 1.5])
    if col_mappable is not None:
        fig.colorbar(
            col_mappable,
            ax=axes.ravel().tolist(),
            cax=cax,
            orientation="vertical",
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


def heatmap_plot_accuracy_metrics_pseudobulk_rd0(
    df, dataset, save_fig=False, FIGURES_PATH=""
):
    metrics = ["fraction_inconclusive", "f1", "precision", "recall"]
    fig_name = f"pseudobulk_{dataset}_accuracy_metrics_heatmap_rd0"
    fontsize = 14

    fig, ax = plt.subplots(2, 2, figsize=(8, 8))
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
                        fontsize=fontsize - 2,
                    )

        ax[i // 2, i % 2].set_title(
            f"{METRIC_LABELS.get(metric, metric)}", fontsize=fontsize
        )

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
            ax[i // 2, i % 2].set_xlabel("# of cells", fontsize=fontsize)
            ax[i // 2, i % 2].set_xticklabels(
                xticklabels, ha="center", fontsize=fontsize
            )
        else:
            ax[i // 2, i % 2].tick_params(axis="x", labelbottom=False)

        # Only show y-label on left column
        if i % 2 == 0:
            ax[i // 2, i % 2].set_ylabel("Tools", fontsize=fontsize)
            ax[i // 2, i % 2].set_yticklabels([TOOL_LABELS.get(label, label) for label in matrix.index], fontsize=fontsize)
        else:
            ax[i // 2, i % 2].tick_params(axis="y", labelleft=False)
    fig.suptitle(
        f"{dataset.replace('_', ' ')} pseudobulks - performance metrics".title(),
        fontsize=fontsize,
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


def heatmap_plot_accuracy_metrics_pseudobulk_rd0_ncells500k(
    df, save_fig=False, FIGURES_PATH=""
):
    metrics = ["fraction_inconclusive", "f1", "precision", "recall"]

    fig_name = f"pseudobulk_accuracy_metrics_heatmap_rd0_ncells500k"
    fontsize = 12

    fig, ax = plt.subplots(2, 2, figsize=(8, 8))
    cmap = plt.get_cmap("Blues")
    cmap.set_bad(color="lightgrey")

    for i, metric in enumerate(metrics):
        matrix = df[(df["read depth"] == 0) & (df["ncells"] == 500000)].pivot_table(
            index=["tool"], columns="dataset", values=metric
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
                        fontsize=fontsize - 2,
                    )

        ax[i // 2, i % 2].set_title(
            f"{METRIC_LABELS.get(metric, metric)}", fontsize=fontsize
        )
        ax[i // 2, i % 2].set_yticks(np.arange(len(matrix.index)))
        
        # Only show x-label on bottom row
        if i // 2 == 1:
            ax[i // 2, i % 2].set_xticks(np.arange(len(matrix.columns)))
            ax[i // 2, i % 2].set_xlabel("Datasets", fontsize=fontsize)
            ax[i // 2, i % 2].set_xticklabels(
                [DATASET_LABELS_SHORT.get(label, label) for label in matrix.columns], ha="right", fontsize=fontsize, rotation=45
            )
        else:
            ax[i // 2, i % 2].tick_params(axis="x", labelbottom=False)

        # Only show y-label on left column
        if i % 2 == 0:
            ax[i // 2, i % 2].set_ylabel("Tools", fontsize=fontsize)
            ax[i // 2, i % 2].set_yticklabels([TOOL_LABELS.get(label, label) for label in matrix.index], fontsize=fontsize)
        else:
            ax[i // 2, i % 2].tick_params(axis="y", labelleft=False)
    fig.suptitle(f"pseudobulks - performance metrics".title(), fontsize=fontsize + 4)

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
    group_colors = plt.get_cmap("viridis")(np.linspace(0, 1, len(n_samples_vals)))

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
                color=group_colors[i],
                ecolor=group_colors[i],
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

            cax = _plot_tool_heatmap(
                ax,
                matrix,
                tool,
                method="imshow",
                symmetric_limit=super_extreme_point,
            )
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


def barplot_matches_nonmatches_na(
    DATA_PATH,
    tools,
    dataset,
    read_depths_mod1,
    read_depths_mod2=None,
    mod1="bulk",
    mod2="single-cell",
    save_fig=False,
    FIGURES_PATH="",
):

    fig_name = f"barplot_matches_{dataset}_{mod1}_vs_{mod2}"

    results_df = count_matches_real_data(
        DATA_PATH=DATA_PATH,
        tools=tools,
        dataset=dataset,
        read_depths_mod1=read_depths_mod1,
        read_depths_mod2=read_depths_mod2,
        mod1=mod1,
        mod2=mod2,
    )

    if read_depths_mod2 is None:
        read_depths_mod2 = read_depths_mod1

    # Stacked barplot of matches, non-matches, and NA counts for each read depth
    n_bars_expected = 1
    n_bars_tool = len(read_depths_mod1)
    width_ratios = [n_bars_expected] + [n_bars_tool] * len(tools)
    total_width = (n_bars_expected + len(tools) * n_bars_tool) * 0.5

    fig, axes = plt.subplots(
        1,
        len(tools) + 1,
        figsize=(total_width, 3.5),
        sharey=True,
        gridspec_kw={"width_ratios": width_ratios},
    )

    # Plot expected results in the first subfigure
    expected_df = results_df[results_df["rd_mod1"] == "expected"].iloc[0]
    axes[0].bar([0], expected_df["matches"], color="orangered", label="Matches")
    axes[0].bar(
        [0],
        expected_df["non_matches"],
        bottom=expected_df["matches"],
        color="lightsalmon",
        label="Non-matches",
    )
    axes[0].bar(
        [0],
        expected_df["NA"],
        bottom=expected_df["matches"] + expected_df["non_matches"],
        color="gray",
        label="NA",
    )
    axes[0].set_xticks([0])
    axes[0].set_xticklabels(["expected"])
    axes[0].set_title("Expected")
    axes[0].set_xlabel("Read depth")
    axes[0].set_ylabel("Count")

    for i, tool in enumerate(tools):
        tool_df = results_df[results_df["tool"] == tool].copy()
        tool_df["rd_label"] = tool_df.apply(
            lambda row: f"{row['rd_mod1']}, {row['rd_mod2']}", axis=1
        )
        tool_df.set_index("rd_label", inplace=True)

        x = np.arange(len(tool_df.index))
        labels = tool_df.index.astype(str)

        axes[i + 1].bar(x, tool_df["matches"], color="orangered", label="Matches")
        axes[i + 1].bar(
            x,
            tool_df["non_matches"],
            bottom=tool_df["matches"],
            color="lightsalmon",
            label="Non-matches",
        )
        axes[i + 1].bar(
            x,
            tool_df["NA"],
            bottom=tool_df["matches"] + tool_df["non_matches"],
            color="gray",
            label="NA",
        )

        axes[i + 1].set_xticks(x)
        axes[i + 1].set_xticklabels(labels)
        axes[i + 1].set_title(tool)
        axes[i + 1].set_xlabel("Read depth")

    fig.suptitle(f"Sample matching results for {dataset} dataset")
    plt.tight_layout()

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


def plot_combined_heatmaps_datasets(
        DATA_PATH,
        FIGURES_PATH,
        save_fig=False,
):  
    fig_name = f"combined_heatmaps_WT_HGSOC_LGG"
    fontsize = 12
    cmap_oranges = plt.get_cmap("Oranges")
    cmap_oranges.set_bad(color="lightgrey")

    # Load and prepare Wilms' tumor heatmap
    mod1_wt = "bulk"
    mod2_wt = "single-nucleus"
    matrix_wt = load_heatmap_data(
        DATA_PATH,
        False,
        "HYSYS",
        "wilms_tumor",
        "null",
        "bulk_0_single-nucleus_0",
        mod1_wt,
        mod2_wt,
    )
    assert matrix_wt is not None, "Matrix for Wilms' tumor is None"
    matrix_wt = _modality_submatrix(matrix_wt, mod1_wt, mod2_wt)
    expected_matches_wt = load_expected_matches_real_data(DATA_PATH, "wilms_tumor")
    try:
        filtered_matrix_wt = order_matrix_by_expected_matches(
            matrix_wt,
            expected_matches_wt,
            mod1_wt,
            mod2_wt,
        )
    except ValueError as exc:
        print(
            "Could not order by expected matches for requested modalities "
            f"({mod1_wt} vs {mod2_wt}): {exc}. Falling back to direct modality filtering."
        )
        filtered_matrix_wt = _subset_matrix_by_modalities(matrix_wt, mod1_wt, mod2_wt)

    # Load and prepare HGSOC heatmap
    mod1_hgsoc = "bulk_dissociated_ribo"
    mod2_hgsoc = "single-cell"
    matrix_hgsoc = load_heatmap_data(
        DATA_PATH,
        False,
        "HYSYS",
        "hgsoc",
        "null",
        "bulk_dissociated_ribo_0_single-cell_0",
        mod1_hgsoc,
        mod2_hgsoc,
    )
    assert matrix_hgsoc is not None, "Matrix for HGSOC is None"
    matrix_hgsoc = _modality_submatrix(matrix_hgsoc, mod1_hgsoc, mod2_hgsoc)
    expected_matches_hgsoc = load_expected_matches_real_data(DATA_PATH, "hgsoc")
    try:
        filtered_matrix_hgsoc = order_matrix_by_expected_matches(
            matrix_hgsoc,
            expected_matches_hgsoc,
            mod1_hgsoc,
            mod2_hgsoc,
        )
    except ValueError as exc:
        print(
            "Could not order by expected matches for requested modalities "
            f"({mod1_hgsoc} vs {mod2_hgsoc}): {exc}. Falling back to direct modality filtering."
        )
        filtered_matrix_hgsoc = _subset_matrix_by_modalities(
            matrix_hgsoc, mod1_hgsoc, mod2_hgsoc
        )

    # Load and prepare low-grade glioma heatmap
    mod1_lgg = "bulk"
    mod2_lgg = "single-cell"
    matrix_lgg = load_heatmap_data(
        DATA_PATH,
        False,
        "CrosscheckFingerprints",
        "low_grade_glioma",
        "null",
        "bulk_0_single-cell_0",
        mod1_lgg,
        mod2_lgg,
    )
    assert matrix_lgg is not None, "Matrix for low-grade glioma is None"
    matrix_lgg = _modality_submatrix(matrix_lgg, mod1_lgg, mod2_lgg)
    expected_matches_lgg = load_expected_matches_real_data(DATA_PATH, "low_grade_glioma")
    try:
        filtered_matrix_lgg = order_matrix_by_expected_matches(
            matrix_lgg,
            expected_matches_lgg,
            mod1_lgg,
            mod2_lgg,
        )
    except ValueError as exc:
        print(
            "Could not order by expected matches for requested modalities "
            f"({mod1_lgg} vs {mod2_lgg}): {exc}. Falling back to direct modality filtering."
        )
        filtered_matrix_lgg = _subset_matrix_by_modalities(matrix_lgg, mod1_lgg, mod2_lgg)
    # Hide expected-but-missing samples in the LGG panel.
    filtered_matrix_lgg = filtered_matrix_lgg.dropna(axis=0, how="all").dropna(
        axis=1, how="all"
    )
    extreme_point = max(
        abs(filtered_matrix_lgg.min().min()), abs(filtered_matrix_lgg.max().max())
    )
    cmap_rd_bu = plt.get_cmap("RdBu")
    cmap_rd_bu.set_bad(color="lightgrey")

    # Keep image aspect ratios and enforce same bottom-row height by scaling bottom column widths.
    # width ratio is proportional to each matrix aspect (n_cols / n_rows).
    hgsoc_aspect = filtered_matrix_hgsoc.shape[1] / max(filtered_matrix_hgsoc.shape[0], 1)
    lgg_aspect = filtered_matrix_lgg.shape[1] / max(filtered_matrix_lgg.shape[0], 1)

    fig = plt.figure(figsize=(10, 12))
    gs = fig.add_gridspec(
        2,
        2,
        height_ratios=[1, 1],
        width_ratios=[hgsoc_aspect, lgg_aspect],
        hspace=-0.1,
        wspace=0.2,
    )

    # Plot Wilms' tumor heatmap
    ax1 = fig.add_subplot(gs[0, :])
    ax1.imshow(filtered_matrix_wt, cmap=cmap_oranges)
    ax1.set_title("Wilms' Tumor", fontsize=fontsize + 2)
    ax1.set_xlabel("Single-nucleus", fontsize=fontsize)
    ax1.set_ylabel("Bulk", fontsize=fontsize)
    ax1.set_xticks(np.arange(len(filtered_matrix_wt.columns)))
    ax1.set_yticks(np.arange(len(filtered_matrix_wt.index)))
    ax1.set_xticklabels([""] * len(filtered_matrix_wt.columns))
    ax1.set_yticklabels([""] * len(filtered_matrix_wt.index))

    # Plot HGSOC heatmap
    ax2 = fig.add_subplot(gs[1, 0])
    ax2.imshow(filtered_matrix_hgsoc, cmap=cmap_oranges)
    ax2.set_title("HGSOC", fontsize=fontsize + 2)
    ax2.set_xlabel("Single-cell", fontsize=fontsize)
    ax2.set_ylabel("rRNA- dissociated bulk", fontsize=fontsize)
    ax2.set_xticks(np.arange(len(filtered_matrix_hgsoc.columns)))
    ax2.set_yticks(np.arange(len(filtered_matrix_hgsoc.index)))
    ax2.set_xticklabels([""] * len(filtered_matrix_hgsoc.columns))
    ax2.set_yticklabels([""] * len(filtered_matrix_hgsoc.index))

    # Plot low-grade glioma heatmap
    ax3 = fig.add_subplot(gs[1, 1])
    ax3.imshow(
        filtered_matrix_lgg,
        cmap=cmap_rd_bu,
        vmin=-extreme_point,
        vmax=extreme_point,
    )
    ax3.set_title("Low-Grade Glioma", fontsize=fontsize + 2)
    ax3.set_xlabel("Single-cell", fontsize=fontsize)
    ax3.set_ylabel("Bulk", fontsize=fontsize)
    ax3.set_xticks(np.arange(len(filtered_matrix_lgg.columns)))
    ax3.set_yticks(np.arange(len(filtered_matrix_lgg.index)))
    ax3.set_xticklabels([""] * len(filtered_matrix_lgg.columns))
    ax3.set_yticklabels([""] * len(filtered_matrix_lgg.index))

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


def row_heatmap_plot_real_data(
    DATA_PATH,
    FIGURES_PATH,
    tools,
    dataset,
    rd_mod1,
    rd_mod2=None,
    mod1="bulk",
    mod2="single-cell",
    hide_missing_samples=False,
    save_fig=False,
):
    """
    Create a row of heatmaps for each tool, showing the similarity matrices for the given dataset and modalities.
    Each heatmap corresponds to a different tool.

    input:
        - DATA_PATH: path to data
        - FIGURES_PATH: path to save figures
        - tools: list of tools to plot
        - dataset: name of dataset
        - rd_mod1: read depth for modality 1
        - rd_mod2: read depth for modality 2 (if None, use rd_mod1)
        - mod1: name of modality 1
        - mod2: name of modality 2
        - save_fig: whether to save the figure

    output:
        - row of heatmaps is shown (or saved if save_fig=True)
    """
    if rd_mod2 is None:
        rd_mod2 = rd_mod1

    fig_name = f"row_heatmaps_{dataset}_{mod1}_rd_{rd_mod1}_vs_{mod2}_rd_{rd_mod2}"
    fontsize = 14

    prepared_matrices: list = []
    width_ratios = []
    for tool in tools:
        matrix = load_heatmap_data(
            DATA_PATH,
            False,
            tool,
            dataset,
            "null",
            f"{mod1}_{rd_mod1}_{mod2}_{rd_mod2}",
            mod1,
            mod2,
        )
        if matrix is None:
            prepared_matrices.append(None)
            width_ratios.append(1)
            continue

        # Filter matrix to only include the requested modalities
        matrix_filtered = _subset_matrix_by_modalities(matrix, mod1, mod2)

        # Order the matrix by expected matches if available
        expected_matches = load_expected_matches_real_data(DATA_PATH, dataset)
        try:
            matrix_filtered = order_matrix_by_expected_matches(
                matrix_filtered,
                expected_matches,
                mod1,
                mod2,
            )
        except ValueError as exc:
            print(
                "Could not order by expected matches for requested modalities "
                f"({mod1} vs {mod2}): {exc}. Falling back to direct modality filtering."
            )

        if matrix_filtered.empty:
            prepared_matrices.append(None)
            width_ratios.append(1)
            continue

        crosscheck_extreme_point = max(
            abs(matrix_filtered.min().min()), abs(matrix_filtered.max().max())
        )
        if hide_missing_samples:
            matrix_filtered = matrix_filtered.dropna(axis=0, how="all").dropna(
                axis=1, how="all"
            )

        prepared_matrices.append(
            {
                "tool": tool,
                "matrix": matrix_filtered,
                "symmetric_limit": crosscheck_extreme_point,
            }
        )
        n_rows = max(matrix_filtered.shape[0], 1)
        n_cols = max(matrix_filtered.shape[1], 1)
        width_ratios.append(n_cols / n_rows)

    fig_height = 4
    fig_width = max(4 * len(tools), fig_height * sum(width_ratios))
    fig, axes = plt.subplots(
        1,
        len(tools),
        figsize=(fig_width, fig_height),
        sharey=True,
        gridspec_kw={"width_ratios": width_ratios},
    )
    if len(tools) == 1:
        axes = [axes]

    fig.subplots_adjust(wspace=0.05)
    fig.suptitle(f"{DATASET_LABELS.get(dataset, dataset)} Similarity Matrices", fontsize=fontsize+2, y=1)

    for i, tool in enumerate(tools):
        ax = axes[i]
        prepared = prepared_matrices[i]
        if prepared is None:
            ax.axis("off")
            ax.text(
                0.5,
                0.5,
                "No data",
                ha="center",
                va="center",
                transform=ax.transAxes,
                fontsize=fontsize,
            )
            continue

        matrix_filtered = prepared["matrix"]
        cax = _plot_tool_heatmap(
            ax,
            matrix_filtered,
            tool,
            method="imshow",
            symmetric_limit=prepared["symmetric_limit"],
        )
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(f"{TOOL_LABELS.get(tool, tool)}", fontsize=fontsize)
        ax.set_xticks(np.arange(len(matrix_filtered.columns)))
        ax.set_yticks(np.arange(len(matrix_filtered.index)))
        ax.set_xticklabels([""] * len(matrix_filtered.columns))
        ax.set_yticklabels([""] * len(matrix_filtered.index))
        ax.set_xlabel(f"{MODALITY_LABELS.get(mod2, mod2)}", fontsize=fontsize)
        if i == 0:
            ax.set_ylabel(f"{MODALITY_LABELS.get(mod1, mod1)}", fontsize=fontsize)

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
