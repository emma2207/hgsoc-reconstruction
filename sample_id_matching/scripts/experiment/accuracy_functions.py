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
    parse_sample_matching_results_crosscheckfingerprints,
    parse_sample_matching_results_hysys,
    parse_sample_matching_results_ngscheckmate,
    parse_heatmap_matrix_vireo,
    parse_heatmap_matrix_hysys,
    true_matches_pseudobulk,
    create_pseudobulk_submatrix_vireo,
    load_sample_matching_results,
    matches_matrix_to_pair_df,
)


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
                            print(
                                f"Error! Do not recognize tool {tool}. Choose one of BAMixChecker, CrosscheckFingerprints, HYSYS, NGSCheckmate, or Vireo."
                            )
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
            DATA_PATH, f"../metadata/sample_matches/sample_matches_{dataset}.csv"
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
    # Make sure all entries are strings
    expected_matches = expected_matches.astype(str)

    return expected_matches


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
                if tool == "CrosscheckFingerprints":
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
                    print(
                        f"Error! Do not recognize tool {tool}. Choose one of CrosscheckFingerprints, HYSYS, NGSCheckmate, or Vireo."
                    )
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
    DATA_PATH, tools, dataset, read_depths, mod1="bulk", mod2="single-cell"
):
    """
    Count the number of matches, non-matches, and NAs for real data sample matching results for each tool and read depth.

    input:
        - DATA_PATH: base path to the data directory
        - tools: list of tools to evaluate (e.g., ["CrosscheckFingerprints", "HYSYS", "NGSCheckmate", "Vireo"])
        - dataset: the dataset name
        - read_depths: list of read depths to evaluate
        - mod1: string to identify the first modality in the sample names (default "bulk")
        - mod2: string to identify the second modality in the sample names (default "single-cell")

    output:
        - results_df: dataframe with the number of matches, non-matches, and NAs for each tool, dataset, and read depth.
    """

    expected_matches = load_expected_matches_real_data(DATA_PATH, dataset)
    mod1_col = [col for col in expected_matches.columns if mod1 in col][0]
    mod2_col = [col for col in expected_matches.columns if mod2 in col][0]
    n_mod1_samples = len(expected_matches[mod1_col].dropna())
    n_mod2_samples = len(expected_matches[mod2_col].dropna())
    n_matches = min(n_mod1_samples, n_mod2_samples)
    n_pairs = n_mod1_samples * n_mod2_samples

    # Start with expected matches
    results = [
        {
            "dataset": dataset,
            "rd": "expected",
            "matches": n_matches,
            "non_matches": n_pairs - n_matches,
            "NA": 0,
        }
    ]
    for tool in tools:
        for rd in read_depths:
            inferred_matches = load_sample_matching_results(
                DATA_PATH=DATA_PATH,
                pseudobulk=False,
                tool=tool,
                dataset=dataset,
                ncells="null",
                rd=rd,
                mod1=mod1,
                mod2=mod2,
            )
            if inferred_matches is None:
                print(f"Warning: No results found for rd={rd}. Skipping.")
                continue
            # Filter to only bulk vs single-cell comparisons
            inferred_matches = inferred_matches.loc[
                [idx for idx in inferred_matches.index if mod1 in idx],
                [col for col in inferred_matches.columns if mod2 in col],
            ]
            # Count matches, non-matches, and NAs
            inferred_matches.shape
            total_pairs = inferred_matches.shape[0] * inferred_matches.shape[1]
            match_counts = (inferred_matches == 1).sum().sum()
            nonmatch_counts = (inferred_matches == 0).sum().sum()
            na_counts = (inferred_matches.isna()).sum().sum()
            if match_counts + nonmatch_counts + na_counts != total_pairs:
                print(f"Warning: Counts do not sum up to total pairs for rd={rd}.")
            results.append(
                {
                    "tool": tool,
                    "dataset": dataset,
                    "rd": rd,
                    "matches": match_counts,
                    "non_matches": nonmatch_counts,
                    "NA": na_counts,
                }
            )
    return pd.DataFrame(results)


def find_mismatches_real_data(
    DATA_PATH, tool, dataset, read_depth, mod1="bulk", mod2="single-cell"
):
    """
    Compare matching sample pairs between expected and inferred matches.

    input:
        - DATA_PATH: base path to the data directory
        - tool: the tool for which to compare expected vs inferred matches (e.g., "CrosscheckFingerprints", "HYSYS", "NGSCheckmate", "Vireo")
        - dataset: the dataset name
        - read_depth: the read depth filter cut-off used for the analysis
        - mod1: string to identify the first modality in the sample names (default "bulk")
        - mod2: string to identify the second modality in the sample names (default "single-cell")

    output:
        - merged_df: dataframe that merges expected matches with inferred matches,
            with an additional column "match_status" indicating whether a pair of
            samples matches is found in the expected matches, inferred matches, or both.
    """

    # Load expected matches
    expected_matches = load_expected_matches_real_data(DATA_PATH, dataset)
    # Load inferred matches
    inferred_matches = load_sample_matching_results(
        DATA_PATH,
        False,
        tool,
        dataset,
        "null",
        read_depth,
        mod1=mod1,
        mod2=mod2,
    )
    if inferred_matches is None:
        print(
            f"No inferred matches found for {tool} on {dataset} at read depth {read_depth}."
        )
        return None, None

    # Convert inferred matches matrix to pairwise dataframe
    matches_df = matches_matrix_to_pair_df(inferred_matches)

    # Split the sample IDs into their components (e.g., "sample1_bulk" -> "sample1", "bulk")
    matches_df[["sample_id_0", "modality_0"]] = matches_df["sample_id_0"].str.rsplit(
        "_", n=1, expand=True
    )
    matches_df[["sample_id_1", "modality_1"]] = matches_df["sample_id_1"].str.rsplit(
        "_", n=1, expand=True
    )

    # Merge with expected matches to determine which inferred matches are correct
    merged_df = expected_matches.merge(
        matches_df,
        right_on=["sample_id_0", "sample_id_1"],
        left_on=[mod1, mod2],
        how="outer",
        indicator="match_status",
    )

    return merged_df
