import numpy as np
import pandas as pd
import os
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    precision_score,
    recall_score,
    f1_score,
)

from analysis_shared_functions import (
    parse_sample_matching_results_bamixchecker,
    parse_sample_matching_results_conpair,
    parse_sample_matching_results_crosscheckfingerprints,
    parse_sample_matching_results_hysys,
    parse_sample_matching_results_ngscheckmate,
    parse_sample_matching_results_ntsm,
    parse_sample_matching_results_omicsprint,
    parse_sample_matching_results_peddy,
    parse_sample_matching_results_somalier,
    parse_sample_matching_results_timeattackgencomp,
    parse_heatmap_matrix_vireo,
    parse_heatmap_matrix_hysys,
    true_matches_pseudobulk,
    create_pseudobulk_submatrix_vireo,
    load_sample_matching_results,
    matches_matrix_to_pair_df,
    _modality_in_label,
)


def _mark_duplicate_value_patterns_inconclusive(inferred_matches):
    """Replace entries in duplicate rows/columns with NaN so they count as inconclusive."""
    duplicate_column_mask = inferred_matches.T.duplicated(keep=False)
    duplicate_row_mask = inferred_matches.duplicated(keep=False)

    if not duplicate_column_mask.any() and not duplicate_row_mask.any():
        return inferred_matches

    inferred_matches = inferred_matches.copy()
    if duplicate_column_mask.any():
        inferred_matches.loc[:, duplicate_column_mask] = np.nan
    if duplicate_row_mask.any():
        inferred_matches.loc[duplicate_row_mask, :] = np.nan

    return inferred_matches


def calculate_accuracy_metrics_pseudobulk(inferred_matches):
    """
    Calculate accuracy metrics (accuracy, balanced accuracy, precision, recall, F1 score) for the inferred sample matches
    compared to the true matches (ground truth) for pseudobulk datasets.

    input:
        - inferred_matches: dataframe with 1 for inferred matches and 0 for non-matches

    output:
        - fraction_inconclusive: the fraction of comparisons that were inconclusive (NaN)
        - accuracy: the overall accuracy of the inferred matches compared to the true matches
        - balanced_accuracy: the balanced accuracy of the inferred matches compared to the true matches
        - precision: the precision of the inferred matches compared to the true matches
        - recall: the recall of the inferred matches compared to the true matches
        - f1: the F1 score of the inferred matches compared to the true matches
    """

    # Calculate fraction of inconclusive (NaN) matches
    total_comparisons = inferred_matches.size
    inconclusive_matches = inferred_matches.isna().sum().sum()
    fraction_inconclusive = inconclusive_matches / total_comparisons

    # Replace NaN with 0 for accuracy calculations (treat inconclusive as no match)
    inferred_matches = inferred_matches.fillna(0).astype(int)

    # Create true matches matrix
    true_matches = true_matches_pseudobulk(inferred_matches)

    # Calculate accuracy metrics by comparing inferred_matches to true_matches
    accuracy = accuracy_score(
        true_matches.values.flatten(), inferred_matches.values.flatten()
    )
    balanced_accuracy = balanced_accuracy_score(
        true_matches.values.flatten(), inferred_matches.values.flatten()
    )
    precision = precision_score(
        true_matches.values.flatten(), inferred_matches.values.flatten()
    )
    recall = recall_score(
        true_matches.values.flatten(), inferred_matches.values.flatten()
    )
    f1 = f1_score(true_matches.values.flatten(), inferred_matches.values.flatten())

    return fraction_inconclusive, accuracy, balanced_accuracy, precision, recall, f1


def loop_accuracy_calculations(
    DATA_PATH, experiment, tools, datasets, ncells_list, read_depths
):
    """
    Loop through all combinations of tools, datasets, number of cells, and read depths to calculate
    accuracy metrics for pseudobulk sample matching.

    input:
        - DATA_PATH: base path to the data directory
        - experiment: the type of experiment (e.g., "missing_samples", "double_samples")
        - tools: list of tools ["BAMixChecker", "CrosscheckFingerprints", "HYSYS", "NGSCheckmate", "Vireo"] to evaluate
        - datasets: list of datasets to evaluate
        - ncells_list: list of number of cells to evaluate
        - read_depths: list of read depths to evaluate

    output:
        - results_df: dataframe with accuracy metrics for each combination of tool, dataset, number of cells, and read depth
    """
    results = []

    for tool in tools:
        for dataset in datasets:
            for ncells in ncells_list:
                for rd in read_depths:
                    try:
                        if tool == "BAMixChecker" and rd == 0:
                            inferred_matches = (
                                parse_sample_matching_results_bamixchecker(
                                    DATA_PATH, True, dataset, ncells
                                )
                            )
                        elif tool == "BAMixChecker" and rd != 0:
                            print(
                                f"Skipping BAMixChecker for read depth {rd} since it does not vary with read depth. "
                            )
                            continue
                        elif tool == "Conpair" and rd == 0:
                            inferred_matches = (
                                parse_sample_matching_results_conpair(
                                    DATA_PATH, True, dataset, ncells
                                )
                            )
                        elif tool == "Conpair" and rd != 0:
                            print(
                                f"Skipping Conpair for read depth {rd} since it does not vary with read depth. "
                            )
                            continue
                        elif tool == "CrosscheckFingerprints":
                            inferred_matches = (
                                parse_sample_matching_results_crosscheckfingerprints(
                                    DATA_PATH, True, dataset, ncells, rd
                                )
                            )
                        elif tool == "HYSYS":
                            inferred_matches = parse_sample_matching_results_hysys(
                                DATA_PATH, True, dataset, ncells, rd
                            )
                        elif tool == "NGSCheckmate":
                            inferred_matches = (
                                parse_sample_matching_results_ngscheckmate(
                                    DATA_PATH, True, dataset, ncells, rd
                                )
                            )
                        elif tool == "ntsm" and rd == 0:
                            inferred_matches = (
                                parse_sample_matching_results_ntsm(
                                    DATA_PATH, True, dataset, ncells,
                                )
                            )
                        elif tool == "ntsm" and rd != 0:
                            print(
                                f"Skipping ntsm for read depth {rd} since it does not vary with read depth. "
                            )
                            continue 
                        elif tool == "OmicsPrint":
                            inferred_matches = (
                                parse_sample_matching_results_omicsprint(
                                    DATA_PATH, True, dataset, ncells, rd
                                )
                            )
                        elif tool == "Peddy":
                            inferred_matches = (
                                parse_sample_matching_results_peddy(
                                    DATA_PATH, True, dataset, ncells, rd,
                                )
                            )
                        elif tool == "Somalier":
                            inferred_matches = (
                                parse_sample_matching_results_somalier(
                                    DATA_PATH, True, dataset, ncells, rd
                                )
                            )
                        elif tool == "TimeAttackGenComp":
                            inferred_matches = (
                                parse_sample_matching_results_timeattackgencomp(
                                    DATA_PATH, True, dataset, ncells, rd
                                )
                            )
                        elif tool == "Vireo":
                            matrix = parse_heatmap_matrix_vireo(
                                DATA_PATH, True, dataset, ncells, rd
                            )
                            _, hysys_matrix = parse_heatmap_matrix_hysys(
                                DATA_PATH, True, dataset, ncells, rd
                            )
                            inferred_matches = create_pseudobulk_submatrix_vireo(
                                matrix, hysys_matrix, experiment
                            )
                        else:
                            print(f"Error! Do not recognize tool {tool}.")
                            print("Choose one of BAMixChecker, Conpair, CrosscheckFingerprints, HYSYS, NGSCheckmate, OmicsPrint, Peddy, Somalier, TimeAttackGenComp, or Vireo.")
                            continue
                    except FileNotFoundError:
                        print(
                            f"Could not find data for {tool}, pseudobulk dataset {dataset}, ncells {ncells}, read depth {rd}. "
                        )
                        continue

                    if experiment == "double_samples":
                        inferred_matches_filtered = inferred_matches.loc[
                            [
                                idx
                                for idx in inferred_matches.index
                                if idx.endswith("_1")
                            ],
                            [
                                col
                                for col in inferred_matches.columns
                                if (col.endswith("_2") or col.endswith("_3"))
                            ],
                        ]
                    elif experiment == "pseudobulk_vs_sc":
                        inferred_matches_filtered = inferred_matches.loc[
                            [
                                idx
                                for idx in inferred_matches.index
                                if idx.endswith("_1")
                            ],
                            [
                                col
                                for col in inferred_matches.columns
                                if col.endswith("_single-cell")
                            ],
                        ]
                    else:
                        inferred_matches_filtered = inferred_matches.loc[
                            [
                                idx
                                for idx in inferred_matches.index
                                if idx.endswith("_1")
                            ],
                            [
                                col
                                for col in inferred_matches.columns
                                if col.endswith("_2")
                            ],
                        ]

                    if tool == "Vireo":
                        inferred_matches_filtered = (
                            _mark_duplicate_value_patterns_inconclusive(
                                inferred_matches_filtered
                            )
                        )


                    if inferred_matches_filtered.size > 0:
                        (
                            fraction_inconclusive,
                            accuracy,
                            balanced_accuracy,
                            precision,
                            recall,
                            f1,
                        ) = calculate_accuracy_metrics_pseudobulk(
                            inferred_matches_filtered
                        )
                    else:
                        fraction_inconclusive = np.nan
                        accuracy = np.nan
                        balanced_accuracy = np.nan
                        precision = np.nan
                        recall = np.nan
                        f1 = np.nan

                    result = {
                        "tool": tool,
                        "pseudobulk": True,
                        "dataset": dataset,
                        "ncells": ncells,
                        "read depth": rd,
                        "fraction_inconclusive": fraction_inconclusive,
                        "accuracy": accuracy,
                        "balanced_accuracy": balanced_accuracy,
                        "precision": precision,
                        "recall": recall,
                        "f1": f1,
                    }
                    results.append(result)

    return pd.DataFrame(results)


def load_expected_matches_real_data(DATA_PATH, dataset):
    """
    Load expected matches for real data from the sample_matches CSV file, clean up the sample names,
    and sort the dataframe by the column with the most non-NaN values.

    input:
        - DATA_PATH: base path to the data directory
        - dataset: name of the dataset

    output:
        - expected_matches: cleaned and sorted dataframe of expected matches for the real data
    """
    expected_matches = pd.read_csv(
        os.path.join(
            DATA_PATH, f"../../metadata/sample_matches/sample_matches_{dataset}.csv"
        )
    )
    expected_matches = expected_matches.dropna(axis=1, how="all").dropna(
        axis=0, how="all"
    )

    expected_matches.columns = expected_matches.columns.str.lower()

    if dataset == "low_grade_glioma":
        expected_matches = expected_matches.drop_duplicates().reset_index(drop=True)
    else:
        # Sort by whichever modality column has the most non-NaN entries (i.e. most samples)
        sort_col = str(expected_matches.count().idxmax())

        # Sort alphabetically by the most-populated modality column,
        # with NaN values appearing at the end.
        expected_matches = (
            expected_matches.sort_values(
                by=sort_col,
                na_position="last",
                kind="mergesort",
            )
            .drop_duplicates()
            .reset_index(drop=True)
        )
    # Make sure all entries are strings while preserving missing values
    expected_matches = expected_matches.astype("string")

    return expected_matches


def _resolve_expected_modality_columns(expected_matches, mod1, mod2):
    """Resolve expected-match columns for requested modalities with graceful fallbacks."""

    def _norm(text):
        return str(text).lower().replace("-", " ").replace("_", " ").strip()

    columns = list(expected_matches.columns)
    norm_cols = {_norm(col): col for col in columns}
    mod1_key = _norm(mod1)
    mod2_key = _norm(mod2)

    mod1_candidates = [orig for norm, orig in norm_cols.items() if mod1_key in norm]
    mod2_candidates = [orig for norm, orig in norm_cols.items() if mod2_key in norm]

    mod1_col = mod1_candidates[0] if len(mod1_candidates) >= 1 else None
    mod2_col = mod2_candidates[0] if len(mod2_candidates) >= 1 else None

    # Fallback: if one modality resolves and table has exactly two columns, pick the other column.
    if mod1_col is not None and mod2_col is None and len(columns) == 2:
        mod2_col = [col for col in columns if col != mod1_col][0]
    if mod2_col is not None and mod1_col is None and len(columns) == 2:
        mod1_col = [col for col in columns if col != mod2_col][0]

    if mod1_col is None or mod2_col is None:
        raise ValueError(
            "Could not resolve expected match columns for requested modalities "
            f"mod1={mod1}, mod2={mod2}. Available columns: {columns}"
        )

    return mod1_col, mod2_col


def accuracy_metrics_averaged_over_iterations(
    DATA_PATH,
    n_iterations,
    dataset,
    ncells,
    rd,
    n_samples_removed,
    experiment="missing_samples",
):
    """
    Loop through all pseudobulk results for a given dataset, number of cells, and read depth,
    and calculate the mean and standard deviation of accuracy metrics for each tool.
    This is used in particular for experiments with missing samples or double samples.

    input:
        - DATA_PATH: base path to the data directory
        - n_iterations: number of iterations to average over
        - dataset: the dataset name
        - ncells: the number of cells used to create the pseudobulk (or "null" for real data)
        - rd: the read depth filter cut-off used for the analysis
        - n_samples_removed: list of number of samples removed (or double samples added) to evaluate
        - experiment: string indicating the experiment type, either "missing_samples" or "double_samples", used for printing messages

    output:
        - final_df: dataframe with mean and standard deviation of accuracy metrics for each tool, dataset,
            number of cells, and read depth.
    """
    final_df = pd.DataFrame()
    if experiment == "missing_samples":
        folder = "remove_samples_"
    elif experiment == "double_samples":
        folder = "double_samples_"
    else:
        print(
            f"Error! Experiment type {experiment} not recognized. Choose either 'missing_samples' or 'double_samples'."
        )
        return final_df
    for n in n_samples_removed:
        for tool in ["CrosscheckFingerprints", "HYSYS", "NGSCheckmate", "Vireo"]:
            accuracy_df = pd.DataFrame()
            for i in range(1, n_iterations + 1):
                # Check if "remove_samples_{n}/it_1/" directory exists
                if not os.path.exists(DATA_PATH + f"{folder}{n}/it_{i}/"):
                    print(
                        f"No results found for {n} samples removed and iteration {i}. Skipping accuracy calculations for this number of samples removed."
                    )
                    continue
                elif tool == "CrosscheckFingerprints":
                    inferred_matches = (
                        parse_sample_matching_results_crosscheckfingerprints(
                            DATA_PATH + f"{folder}{n}/it_{i}/",
                            True,
                            dataset,
                            ncells,
                            rd,
                        )
                    )
                elif tool == "HYSYS":
                    inferred_matches = parse_sample_matching_results_hysys(
                        DATA_PATH + f"{folder}{n}/it_{i}/", True, dataset, ncells, rd
                    )
                elif tool == "NGSCheckmate":
                    inferred_matches = parse_sample_matching_results_ngscheckmate(
                        DATA_PATH + f"{folder}{n}/it_{i}/", True, dataset, ncells, rd
                    )
                elif tool == "Vireo":
                    matrix = parse_heatmap_matrix_vireo(
                        DATA_PATH + f"{folder}{n}/it_{i}/", True, dataset, ncells, rd
                    )
                    _, hysys_matrix = parse_heatmap_matrix_hysys(
                        DATA_PATH + f"{folder}{n}/it_{i}/", True, dataset, ncells, rd
                    )
                    n_samples_double = n if experiment == "double_samples" else 0
                    inferred_matches = create_pseudobulk_submatrix_vireo(
                        matrix, hysys_matrix, experiment, n_samples_double
                    )
                else:
                    print(f"Error! Do not recognize tool {tool}.")
                    print("Choose one of CrosscheckFingerprints, HYSYS, NGSCheckmate, or Vireo.")
                    continue

                if experiment == "double_samples":
                    inferred_matches_filtered = inferred_matches.loc[
                        [idx for idx in inferred_matches.index if idx.endswith("_1")],
                        [
                            col
                            for col in inferred_matches.columns
                            if (col.endswith("_2") or col.endswith("_3"))
                        ],
                    ]
                else:
                    inferred_matches_filtered = inferred_matches.loc[
                        [idx for idx in inferred_matches.index if idx.endswith("_1")],
                        [col for col in inferred_matches.columns if col.endswith("_2")],
                    ]
                if tool == "Vireo":
                    inferred_matches_filtered = (
                        _mark_duplicate_value_patterns_inconclusive(
                            inferred_matches_filtered
                        )
                    )
                (
                    fraction_inconclusive,
                    accuracy,
                    balanced_accuracy,
                    precision,
                    recall,
                    f1,
                ) = calculate_accuracy_metrics_pseudobulk(inferred_matches_filtered)

                accuracy_df = pd.concat(
                    [
                        accuracy_df,
                        pd.DataFrame(
                            {
                                "iteration": [i],
                                "fraction_inconclusive": [fraction_inconclusive],
                                "accuracy": [accuracy],
                                "balanced_accuracy": [balanced_accuracy],
                                "precision": [precision],
                                "recall": [recall],
                                "f1": [f1],
                            }
                        ),
                    ]
                )

            final_df = pd.concat(
                [
                    final_df,
                    pd.DataFrame(
                        {
                            "tool": [tool],
                            "dataset": [dataset],
                            "ncells": [ncells],
                            "rd": [rd],
                            "n_samples_removed": [n],
                            "av_fraction_inconclusive": accuracy_df[
                                "fraction_inconclusive"
                            ].mean(),
                            "av_accuracy": accuracy_df["accuracy"].mean(),
                            "av_balanced_accuracy": accuracy_df[
                                "balanced_accuracy"
                            ].mean(),
                            "av_precision": accuracy_df["precision"].mean(),
                            "av_recall": accuracy_df["recall"].mean(),
                            "av_f1": accuracy_df["f1"].mean(),
                            "sd_fraction_inconclusive": accuracy_df[
                                "fraction_inconclusive"
                            ].std(),
                            "sd_accuracy": accuracy_df["accuracy"].std(),
                            "sd_balanced_accuracy": accuracy_df[
                                "balanced_accuracy"
                            ].std(),
                            "sd_precision": accuracy_df["precision"].std(),
                            "sd_recall": accuracy_df["recall"].std(),
                            "sd_f1": accuracy_df["f1"].std(),
                        }
                    ),
                ]
            )
    return final_df.reset_index(drop=True)


def accuracy_metrics_hysys_heatmap_vireo_algorithm(
    DATA_PATH,
    dataset,
    ncells,
    rd,
):
    """
    Calculate accuracy metrics for HYSYS's heatmap results using the same algorithm Vireo uses to infer sample matches
    (minimize diagonal of the matrix).

    input:
        - DATA_PATH: base path to the data directory
        - dataset: the dataset name
        - ncells: the number of cells used to create the pseudobulk (or "null" for real data)
        - rd: the read depth filter cut-off used for the analysis

    """

    _, hysys_matrix = parse_heatmap_matrix_hysys(DATA_PATH, True, dataset, ncells, rd)
    # Filter down to samples ending in _1 in rows and samples ending in _2 in columns
    inferred_matches = create_pseudobulk_submatrix_vireo(
        hysys_matrix, hysys_matrix, "Basic"
    )

    (
        fraction_inconclusive,
        accuracy,
        balanced_accuracy,
        precision,
        recall,
        f1,
    ) = calculate_accuracy_metrics_pseudobulk(inferred_matches)

    result = {
        "tool": "HYSYS matrix, Vireo algorithm",
        "pseudobulk": True,
        "dataset": dataset,
        "ncells": ncells,
        "read depth": rd,
        "fraction_inconclusive": fraction_inconclusive,
        "accuracy": accuracy,
        "balanced_accuracy": balanced_accuracy,
        "precision": precision,
        "recall": recall,
        "f1": f1,
    }

    return pd.DataFrame(result, index=[0])


def count_matches_real_data(
    DATA_PATH,
    tools,
    dataset,
    read_depths_mod1,
    read_depths_mod2=None,
    mod1="bulk",
    mod2="single-cell",
):
    """
    Count the number of matches, non-matches, and NAs for real data sample matching
    results for each tool and read-depth pair.

    input:
        - DATA_PATH: base path to the data directory
        - tools: list of tools to evaluate (e.g., ["CrosscheckFingerprints", "HYSYS", "NGSCheckmate", "Vireo"])
        - dataset: the dataset name
        - read_depths_mod1: list of read depths for modality 1. If read_depths_mod2 is
            None, these values are reused for modality 2 (backward compatible behavior).
        - read_depths_mod2: optional list of read depths for modality 2.
        - mod1: string to identify the first modality in the sample names (default "bulk")
        - mod2: string to identify the second modality in the sample names (default "single-cell")

    output:
        - results_df: dataframe with the number of matches, non-matches, and NAs for
            each tool, dataset, and read-depth pair.
    """

    expected_matches = load_expected_matches_real_data(DATA_PATH, dataset)
    mod1_col, mod2_col = _resolve_expected_modality_columns(
        expected_matches, mod1, mod2
    )
    n_mod1_samples = len(expected_matches[mod1_col].dropna())
    n_mod2_samples = len(expected_matches[mod2_col].dropna())
    n_matches = min(n_mod1_samples, n_mod2_samples)
    n_pairs = n_mod1_samples * n_mod2_samples

    # Start with expected matches
    results = [
        {
            "dataset": dataset,
            "rd_mod1": "expected",
            "rd_mod2": "expected",
            "matches": n_matches,
            "non_matches": n_pairs - n_matches,
            "NA": 0,
        }
    ]

    if read_depths_mod2 is None:
        read_depths_mod2 = read_depths_mod1
    if len(read_depths_mod1) != len(read_depths_mod2):
        raise ValueError(
            "read_depths_mod1 and read_depths_mod2 must have the same length."
        )

    for tool in tools:
        for rd1, rd2 in zip(read_depths_mod1, read_depths_mod2):
            inferred_matches = load_sample_matching_results(
                DATA_PATH=DATA_PATH,
                pseudobulk=False,
                tool=tool,
                dataset=dataset,
                ncells="null",
                rd1=rd1,
                rd2=rd2,
                mod1=mod1,
                mod2=mod2,
            )
            if inferred_matches is None:
                print(f"Warning: No results found for rd1={rd1}, rd2={rd2}. Skipping.")
                continue
            # Filter to only bulk vs single-cell comparisons
            inferred_matches = inferred_matches.loc[
                [
                    idx
                    for idx in inferred_matches.index
                    if _modality_in_label(mod1, idx)
                ],
                [
                    col
                    for col in inferred_matches.columns
                    if _modality_in_label(mod2, col)
                ],
            ]
            # Count matches, non-matches, and NAs
            total_pairs = inferred_matches.shape[0] * inferred_matches.shape[1]
            match_counts = (inferred_matches == 1).sum().sum()
            nonmatch_counts = (inferred_matches == 0).sum().sum()
            na_counts = (inferred_matches.isna()).sum().sum()
            if match_counts + nonmatch_counts + na_counts != total_pairs:
                print(
                    f"Warning: Counts do not sum up to total pairs for rd1={rd1}, rd2={rd2}."
                )
            results.append(
                {
                    "tool": tool,
                    "dataset": dataset,
                    "rd_mod1": rd1,
                    "rd_mod2": rd2,
                    "matches": match_counts,
                    "non_matches": nonmatch_counts,
                    "NA": na_counts,
                }
            )
    return pd.DataFrame(results)


def find_mismatches_real_data(
    DATA_PATH, tool, dataset, rd1, rd2, mod1="bulk", mod2="single-cell"
):
    """
    Compare matching sample pairs between expected and inferred matches.

    input:
        - DATA_PATH: base path to the data directory
        - tool: the tool for which to compare expected vs inferred matches (e.g., "CrosscheckFingerprints", "HYSYS", "NGSCheckmate", "Vireo")
        - dataset: the dataset name
        - rd1: read depth filter cut-off for modality 1
        - rd2: read depth filter cut-off for modality 2
        - mod1: string to identify the first modality in the sample names (default "bulk")
        - mod2: string to identify the second modality in the sample names (default "single-cell")

    output:
        - merged_df: dataframe that merges expected matches with inferred matches,
            with an additional column "match_status" indicating whether a pair of
            samples matches is found in the expected matches, inferred matches, or both.
    """

    # Load expected matches
    expected_matches = load_expected_matches_real_data(DATA_PATH, dataset)
    mod1_col, mod2_col = _resolve_expected_modality_columns(
        expected_matches, mod1, mod2
    )
    # Load inferred matches
    inferred_matches = load_sample_matching_results(
        DATA_PATH,
        False,
        tool,
        dataset,
        "null",
        rd1,
        rd2,
        mod1=mod1,
        mod2=mod2,
    )
    if inferred_matches is None:
        print(
            f"No inferred matches found for {tool} on {dataset} at rd1={rd1}, rd2={rd2}."
        )
        return None

    # Convert inferred matches matrix to pairwise dataframe
    matches_df = matches_matrix_to_pair_df(inferred_matches)

    # Standardize pair column names across tools (some matrices keep custom index/column names).
    if (
        "sample_id_0" not in matches_df.columns
        or "sample_id_1" not in matches_df.columns
    ):
        if len(matches_df.columns) >= 2:
            matches_df = matches_df.rename(
                columns={
                    matches_df.columns[0]: "sample_id_0",
                    matches_df.columns[1]: "sample_id_1",
                }
            )
        else:
            raise ValueError(
                "Could not identify sample-pair columns in inferred matches dataframe."
            )

    # Split the sample IDs into their components (e.g., "sample1_bulk" -> "sample1", "bulk")
    # Use temporary variables to avoid pandas multi-column assignment issues
    # Handle cases where modality suffix may not be present (rsplit returns 1 column instead of 2)
    split_0 = matches_df["sample_id_0"].str.rsplit("_", n=1, expand=True)
    if split_0.shape[1] == 2:
        matches_df["modality_0"] = split_0.iloc[:, 1]
        matches_df["sample_id_0"] = split_0.iloc[:, 0]
    else:
        # No underscore found - sample_id stays as-is, modality is None
        matches_df["modality_0"] = None

    split_1 = matches_df["sample_id_1"].str.rsplit("_", n=1, expand=True)
    if split_1.shape[1] == 2:
        matches_df["modality_1"] = split_1.iloc[:, 1]
        matches_df["sample_id_1"] = split_1.iloc[:, 0]
    else:
        # No underscore found - sample_id stays as-is, modality is None
        matches_df["modality_1"] = None

    # Merge with expected matches to determine which inferred matches are correct
    merged_df = expected_matches.merge(
        matches_df,
        right_on=["sample_id_0", "sample_id_1"],
        left_on=[mod1_col, mod2_col],
        how="outer",
        indicator="match_status",
    )

    return merged_df
