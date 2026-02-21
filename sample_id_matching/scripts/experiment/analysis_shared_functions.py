import os
import re
import numpy as np
import pandas as pd

from sklearn.metrics import (
    recall_score,
    precision_score,
    f1_score,
    balanced_accuracy_score,
    accuracy_score,
)
from scipy.optimize import linear_sum_assignment


dataset_regex_dict = {
    "hgsoc": r"[a-zA-Z0-9]{32}_",
    "low_grade_glioma": r"^GSM[0-9]{7}_",
}


####################################################
### Functions to convert long dataframes to matrices
####################################################
def long_df_to_matrix_bamixchecker(df, metric="Concordance Rate"):
    # Use pivot_table for much faster matrix creation
    matrix = df.pivot_table(
        index=0,
        columns="Sample",
        values=metric,
        aggfunc="first",  # Use first value if duplicates exist
    )

    # Get all unique samples to create a square matrix
    all_samples = sorted(set(df.index) | set(df["Sample"]))
    matrix = matrix.reindex(index=all_samples, columns=all_samples)

    # Make matrix symmetric by filling NaN values
    matrix = matrix.combine_first(matrix.T)

    # Order rows and columns
    matrix = matrix[matrix.columns.sort_values()]
    matrix = matrix.reindex(sorted(matrix.index))

    return matrix


def long_df_to_matrix_crosscheckfingerprints(df, metric="LOD_SCORE"):
    # Use pivot_table for much faster matrix creation
    matrix = df.pivot_table(
        index="LEFT_SAMPLE",
        columns="RIGHT_SAMPLE",
        values=metric,
        aggfunc="first",  # Use first value if duplicates exist
    )

    # Get all unique samples to create a square matrix
    all_samples = sorted(set(df["LEFT_SAMPLE"]) | set(df["RIGHT_SAMPLE"]))
    matrix = matrix.reindex(index=all_samples, columns=all_samples)

    # Make matrix symmetric by filling NaN values
    matrix = matrix.combine_first(matrix.T)

    # Order rows and columns by bulk modality and sample ID
    matrix = matrix.sort_index(key=lambda x: x.str.split("_").str[0]).sort_index(
        key=lambda x: x.str.split("_", n=1).str[1]
    )
    matrix = matrix.reindex(
        sorted(
            sorted(matrix.columns, key=lambda x: x.split("_", 1)[0]),
            key=lambda x: x.split("_", 1)[1],
        )
    )
    matrix = matrix[matrix.index]

    return matrix


def long_df_to_matrix_hysys(df):
    # Use pivot_table for much faster matrix creation
    matrix = df.pivot_table(
        index=0,
        columns="Sample",
        values="Concordance",
        aggfunc="first",  # Use first value if duplicates exist
    )

    # Get all unique samples to create a square matrix
    all_samples = sorted(set(df.index) | set(df["Sample"]))
    matrix = matrix.reindex(index=all_samples, columns=all_samples)

    # Make matrix symmetric by filling NaN values
    matrix = matrix.combine_first(matrix.T)

    # Order rows and columns by bulk modality and sample ID
    matrix = matrix.sort_index(key=lambda x: x.str.split("_").str[0]).sort_index(
        key=lambda x: x.str.split("_", n=1).str[1]
    )
    matrix = matrix.reindex(
        sorted(
            sorted(matrix.columns, key=lambda x: x.split("_", 1)[0]),
            key=lambda x: x.split("_", 1)[1],
        )
    )
    matrix = matrix[matrix.index]

    return matrix


def long_df_to_matrix_ngscheckmate(df, metric="Correlation"):
    # Use pivot_table for much faster matrix creation
    matrix = df.pivot_table(
        index=0,
        columns="Sample",
        values=metric,
        aggfunc="first",  # Use first value if duplicates exist
    )

    # Convert correlation values to percentage (divide by 100)
    if metric == "Correlation":
        matrix = matrix / 100

    # Get all unique samples to create a square matrix
    all_samples = sorted(set(df.index) | set(df["Sample"]))
    matrix = matrix.reindex(index=all_samples, columns=all_samples)

    # Make matrix symmetric by filling NaN values
    matrix = matrix.combine_first(matrix.T)

    # Order rows and columns by bulk modality and sample ID
    matrix = matrix.sort_index(key=lambda x: x.str.split("_").str[0]).sort_index(
        key=lambda x: x.str.split("_", n=1).str[1]
    )
    matrix = matrix.reindex(
        sorted(
            sorted(matrix.columns, key=lambda x: x.split("_", 1)[0]),
            key=lambda x: x.split("_", 1)[1],
        )
    )
    matrix = matrix[matrix.index]

    return matrix


#####################################################
### Functions to parse tool outputs and create heatmap matrices
#####################################################
def parse_heatmap_matrix_bamixchecker(DATA_PATH, pseudobulk, dataset, ncells):
    if pseudobulk:
        df = pd.read_csv(
            DATA_PATH
            + f"2b_bamixchecker/{dataset}/pseudobulk/ncells_{ncells}/BAMixChecker/Total_result.txt",
            sep="\t",
            index_col=0,
            header=None,
        )
    else:
        df = pd.read_csv(
            DATA_PATH
            + f"2b_bamixchecker/{dataset}/real_data/ncells_null/BAMixChecker/Total_result.txt",
            sep="\t",
            index_col=0,
            header=None,
        )

    df.columns = ["Sample", "Concordance Rate", "Conclusion"]

    regex_exp = dataset_regex_dict.get(dataset, "")
    rename_dict = {
        x: re.sub(regex_exp, "", x.replace("_filtered_rg.bam", ""))
        for x in list(set(df.index) | set(df["Sample"]))
    }
    df.index = df.index.map(rename_dict)
    df["Sample"] = df["Sample"].map(rename_dict)

    matrix = long_df_to_matrix_bamixchecker(df)

    return df, matrix


def parse_heatmap_matrix_crosscheckfingerprints(
    DATA_PATH, pseudobulk, dataset, ncells, rd
):
    # Read CrosscheckFingerprints output
    if pseudobulk:
        df = pd.read_csv(
            DATA_PATH
            + f"2a_fingerprints/{dataset}/pseudobulk/ncells_{ncells}/read_depth_{rd}/crosscheck_metrics.txt",
            sep="\t",
            skiprows=6,
            header=0,
        )
    else:
        df = pd.read_csv(
            DATA_PATH
            + f"2a_fingerprints/{dataset}/real_data/ncells_null/read_depth_{rd}/crosscheck_metrics.txt",
            sep="\t",
            skiprows=6,
            header=0,
        )
    df = df[["LEFT_SAMPLE", "RIGHT_SAMPLE", "LOD_SCORE", "RESULT"]]

    if pseudobulk:
        regex_exp = dataset_regex_dict.get(dataset, "")
    else:
        regex_exp = ""

    rename_dict = {
        x: re.sub(regex_exp, "", x.replace(".bam", "").replace("ds.", ""))
        for x in list(set(df["LEFT_SAMPLE"]) | set(df["RIGHT_SAMPLE"]))
    }
    df["LEFT_SAMPLE"] = df["LEFT_SAMPLE"].map(rename_dict)
    df["RIGHT_SAMPLE"] = df["RIGHT_SAMPLE"].map(rename_dict)

    # Create matrix for heatmap
    matrix = long_df_to_matrix_crosscheckfingerprints(df, metric="LOD_SCORE")

    return df, matrix


def parse_heatmap_matrix_hysys(DATA_PATH, pseudobulk, dataset, ncells, rd):
    if pseudobulk:
        df = pd.read_csv(
            DATA_PATH
            + f"2a_hysys/{dataset}/pseudobulk/ncells_{ncells}/read_depth_{rd}/concordance_output.txt",
            sep="\t",
            index_col=0,
            header=None,
        )
    else:
        df = pd.read_csv(
            DATA_PATH
            + f"2a_hysys/{dataset}/real_data/ncells_null/read_depth_{rd}/concordance_output.txt",
            sep="\t",
            index_col=0,
            header=None,
        )

    regex_exp = dataset_regex_dict.get(dataset, "")

    rename_dict = {
        x: re.sub(
            regex_exp,
            "",
            x.replace(f"{dataset}/", "")
            .replace("_individual_variants.snps", "")
            .replace("pseudobulk/", "")
            .replace("real_data/", "")
            .replace(f"ncells_{ncells}/", "")
            .replace(f"read_depth_{rd}/", ""),
        )
        for x in list(set(df.index) | set(df[1]))
    }
    df.index = df.index.map(rename_dict)
    df[1] = df[1].map(rename_dict)
    df.columns = ["Sample", "Concordance", "Size"]

    matrix = long_df_to_matrix_hysys(df)

    return df, matrix


def parse_heatmap_matrix_ngscheckmate(DATA_PATH, pseudobulk, dataset, ncells, rd):
    # Read NGSCheckmate output
    if pseudobulk:
        df = pd.read_csv(
            DATA_PATH
            + f"2a_ngscheckmate/{dataset}/pseudobulk/ncells_{ncells}/read_depth_{rd}/output_all.txt",
            sep="\t",
            index_col=0,
            header=None,
        )
    else:
        df = pd.read_csv(
            DATA_PATH
            + f"2a_ngscheckmate/{dataset}/real_data/ncells_{ncells}/read_depth_{rd}/output_all.txt",
            sep="\t",
            index_col=0,
            header=None,
        )

    regex_exp = dataset_regex_dict.get(dataset, "")

    # Clean up sample names
    rename_dict = {
        x: re.sub(
            regex_exp,
            "",
            x.replace(".vcf", "")
            .replace("ds.", "")
            .replace("_individual_variants", ""),
        )
        for x in list(set(df.index) | set(df[2]))
    }
    df.index = df.index.map(rename_dict)
    df[2] = df[2].map(rename_dict)
    # Rename columns
    df.columns = ["Matched", "Sample", "Binary", "Correlation"]

    # Create heatmap matrix
    matrix = long_df_to_matrix_ngscheckmate(df)

    return df, matrix


def parse_heatmap_matrix_vireo(DATA_PATH, pseudobulk, dataset, ncells, rd):
    if pseudobulk:
        matrix = pd.read_csv(
            DATA_PATH
            + f"2a_vireo/{dataset}/pseudobulk/ncells_{ncells}/read_depth_{rd}/similarity_matrix.csv",
            index_col=0,
        )
    else:
        matrix = pd.read_csv(
            DATA_PATH
            + f"2a_vireo/{dataset}/real_data/ncells_{ncells}/read_depth_{rd}/similarity_matrix.csv",
            index_col=0,
        )

    regex_exp = dataset_regex_dict.get(dataset, "")

    matrix.columns = [
        re.sub(regex_exp, "", col.replace(".bam", "").replace("ds.", ""))
        for col in matrix.columns
    ]
    matrix.index = pd.Index(
        [
            re.sub(regex_exp, "", idx.replace(".bam", "").replace("ds.", ""))
            for idx in matrix.index
        ]
    )

    # Order rows and columns
    ordered_columns = sorted(matrix.columns)
    ordered_index = sorted(matrix.index)
    matrix = matrix.loc[ordered_index, ordered_columns]

    return matrix


#######################################################
### Functions to parse tool outputs and create sample match dataframes
#######################################################
def parse_sample_matching_results_bamixchecker(DATA_PATH, pseudobulk, dataset, ncells):

    # Create a mapping for BAMixChecker results
    # 1 = match, 0 = no match
    # Note BAMixChecker does not distinguish between "no match" and "inconclusive"
    result_mapping = {
        "Unmatched": 0,
        "Matched": 1,
    }
    df, _ = parse_heatmap_matrix_bamixchecker(DATA_PATH, pseudobulk, dataset, ncells)
    sample_matches = long_df_to_matrix_bamixchecker(df, metric="Conclusion").replace(
        result_mapping
    )
    return sample_matches


def parse_sample_matching_results_crosscheckfingerprints(
    DATA_PATH, pseudobulk, dataset, ncells, rd
):

    # Create a mapping for CrosscheckFingerprints results
    # 1 = match, 0 = no match, nan = inconclusive
    result_mapping = {
        "EXPECTED_MATCH": 1,
        "UNEXPECTED_MATCH": 1,
        "EXPECTED_MISMATCH": 0,
        "UNEXPECTED_MISMATCH": 0,
        "INCONCLUSIVE": np.nan,
    }
    df, _ = parse_heatmap_matrix_crosscheckfingerprints(
        DATA_PATH, pseudobulk, dataset, ncells, rd
    )
    sample_matches = long_df_to_matrix_crosscheckfingerprints(df, metric="RESULT")
    sample_matches = sample_matches.replace(result_mapping)

    return sample_matches


def parse_sample_matching_results_hysys(DATA_PATH, pseudobulk, dataset, ncells, rd):
    """
    Parse HYSYS model_results.txt file to categorize sample relationships.

    Args:
        file_path (str): Path to the model_results.txt file
        dataset (str): Name of the dataset

    Returns:
        pd.DataFrame: DataFrame with MultiIndex (sample1, sample2) and 'match' column
                     Values: 1 (matching), 0 (not matching), NaN (inconclusive)
    """
    if pseudobulk:
        file_path = os.path.join(
            DATA_PATH,
            f"2a_hysys/{dataset}/pseudobulk/ncells_{ncells}/read_depth_{rd}/model_results.txt",
        )
    else:
        file_path = os.path.join(
            DATA_PATH,
            f"2a_hysys/{dataset}/real_data/ncells_null/read_depth_{rd}/model_results.txt",
        )

    with open(file_path, "r") as file:
        content = file.read()

    # Initialize collections
    matching_pairs = set()
    inconclusive_pairs = set()
    all_samples = set()

    # Clean sample names function
    def clean_sample_name(sample_path):
        # Extract filename from path
        sample_name = sample_path.split("/")[-1]
        # Remove file extension
        sample_name = sample_name.replace(".snps", "")
        # Apply same regex cleanup as other tools
        regex_exp = dataset_regex_dict.get(dataset, "")
        sample_name = re.sub(regex_exp, "", sample_name)
        sample_name = sample_name.replace("_individual_variants", "").replace("ds.", "")
        return sample_name

    # Find all "reject concordance value" lines (inconclusive)
    reject_pattern = r"reject concordance value.*?\(([^,]+),\s*([^)]+)\)"
    reject_matches = re.findall(reject_pattern, content)

    for sample1_path, sample2_path in reject_matches:
        sample1 = clean_sample_name(sample1_path.strip())
        sample2 = clean_sample_name(sample2_path.strip())
        all_samples.update([sample1, sample2])
        if sample1 != sample2:  # Don't include self-comparisons
            pair = tuple(sorted([sample1, sample2]))  # Sort to avoid duplicates
            inconclusive_pairs.add(pair)

    # Find all "sample X is related to: [Y]" lines (matching)
    related_pattern = r"sample ([^\s]+) is related to: \[([^\]]+)\]"
    related_matches = re.findall(related_pattern, content)

    for sample1_path, related_samples_str in related_matches:
        sample1 = clean_sample_name(sample1_path.strip())
        all_samples.add(sample1)

        # Parse the list of related samples
        related_samples = [
            s.strip().strip("'\"") for s in related_samples_str.split(",")
        ]

        for sample2_path in related_samples:
            sample2 = clean_sample_name(sample2_path.strip())
            all_samples.add(sample2)
            if sample1 != sample2:  # Don't include self-relationships
                pair = tuple(sorted([sample1, sample2]))  # Sort to avoid duplicates
                matching_pairs.add(pair)

    # Create all possible pairs from all samples
    all_pairs = [(s1, s2) for s1 in all_samples for s2 in all_samples]

    # Create DataFrame with all pairs
    pair_data = []
    for pair in all_pairs:
        if pair in matching_pairs:
            match_status = 1  # matching
        elif pair in inconclusive_pairs:
            match_status = float("nan")  # inconclusive
        else:
            match_status = 0  # not matching

        pair_data.append(
            {"sample1": pair[0], "sample2": pair[1], "match": match_status}
        )

    # Create matrix with sample1 as index and sample2 as columns
    df = pd.DataFrame(pair_data)
    matrix = df.pivot_table(index="sample1", columns="sample2", values="match")

    return matrix


def parse_sample_matching_results_ngscheckmate(
    DATA_PATH, pseudobulk, dataset, ncells, rd
):

    result_mapping = {
        "unmatched": 0,
        "matched": 1,
    }
    df, _ = parse_heatmap_matrix_ngscheckmate(
        DATA_PATH, pseudobulk, dataset, ncells, rd
    )
    sample_matches = long_df_to_matrix_ngscheckmate(df, "Matched").replace(
        result_mapping
    )

    return sample_matches


def parse_sample_matching_results_vireo(DATA_PATH, pseudobulk, dataset, ncells, rd):
    if pseudobulk:
        sample_matches = pd.read_csv(
            DATA_PATH
            + f"2a_vireo/{dataset}/pseudobulk/ncells_{ncells}/read_depth_{rd}/matched_samples.csv",
            index_col=0,
        ).reset_index()
    else:
        sample_matches = pd.read_csv(
            DATA_PATH
            + f"2a_vireo/{dataset}/real_data/ncells_{ncells}/read_depth_{rd}/matched_samples.csv",
            index_col=0,
        ).reset_index()

    regex_exp = dataset_regex_dict.get(dataset, "")
    for col in sample_matches.columns:
        sample_matches[col] = [
            re.sub(regex_exp, "", sample.replace(".bam", "").replace("ds.", ""))
            for sample in sample_matches[col]
        ]
    sample_matches = sample_matches.sort_values(by=["sample_id_0"]).reset_index(
        drop=True
    )

    # Vireo will not include all samples in the results if we're comparing two sets of samples of unequal length
    # We have to get a list of all samples from another tool's results
    hysys_df, _ = parse_heatmap_matrix_hysys(DATA_PATH, pseudobulk, dataset, ncells, rd)
    all_samples = sorted(set(hysys_df.index))

    # Create all possible pairs of samples
    all_pairs = [(s1, s2) for s1 in all_samples for s2 in all_samples]
    all_pairs_df = pd.DataFrame(all_pairs, columns=["sample_id_0", "sample_id_1"])

    # Add a "match" column: 1 if pair is in sample_matches, 0 otherwise
    sample_matches_set = set(
        zip(sample_matches["sample_id_0"], sample_matches["sample_id_1"])
    )
    all_pairs_df["match"] = all_pairs_df.apply(
        lambda row: (
            1 if (row["sample_id_0"], row["sample_id_1"]) in sample_matches_set else 0
        ),
        axis=1,
    )

    return all_pairs_df.set_index(["sample_id_0", "sample_id_1"])


######################################################
### Functions to calculate accuracy metrics
######################################################
def create_pseudobulk_submatrix_vireo(matrix):
    # Filter Vireo matrix to only samples ending in _1 in rows and samples ending in _2 in columns (or _2 vs _3 or _1 vs _3)
    pseudobulk_1 = "_1"
    pseudobulk_2 = "_2"
    matrix_filtered = matrix.loc[
        [idx for idx in matrix.index if idx.endswith(pseudobulk_1)],
        [col for col in matrix.columns if col.endswith(pseudobulk_2)],
    ]

    # Minimize diagonal of matrix_filtered by moving rows and columns around
    # This is the algorith Vireo uses to match samples
    idx0, idx1 = linear_sum_assignment(matrix_filtered.values)
    matrix_reordered = pd.DataFrame(
        matrix_filtered.iloc[idx0, idx1],
        index=matrix_filtered.index[idx0],
        columns=matrix_filtered.columns[idx1],
    )
    # The new order of samples in the index and columns reveals the matching of samples between the two pseudobulk sets.
    # Create a sample-matching matrix with 1 on the diagonal and 0 elsewhere, where the order of rows and columns is the same as in matrix_reordered.
    # This matrix has the inferred sample matches
    inferred_matches = pd.DataFrame(
        np.identity(len(matrix_reordered.index)),
        index=pd.Index(matrix_reordered.index),
        columns=pd.Index(matrix_reordered.columns),
    )
    # Reorder to the original order of samples
    ordered_columns = sorted(inferred_matches.columns)
    ordered_index = sorted(inferred_matches.index)
    inferred_matches = inferred_matches.loc[ordered_index, ordered_columns]

    return inferred_matches


def true_matches_pseudobulk(matrix):
    # Make sure the columns are ordered
    matrix = matrix.loc[sorted(matrix.index), sorted(matrix.columns)]

    # Create matrix with true matches
    true_matches = pd.DataFrame(
        np.identity(len(matrix)),
        index=pd.Index(matrix.index),
        columns=pd.Index(matrix.columns),
    )

    return true_matches


def calculate_accuracy_metrics_pseudobulk(inferred_matches):

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


def loop_accuracy_calculations(DATA_PATH, tools, datasets, ncells_list, read_depths):
    results = []

    for tool in tools:
        for dataset in datasets:
            for ncells in ncells_list:
                for rd in read_depths:
                    try:
                        if (tool == "BAMixChecker" and rd == 1):
                            inferred_matches = (
                                parse_sample_matching_results_bamixchecker(
                                    DATA_PATH, True, dataset, ncells
                                )
                            )
                        elif tool == "BAMixChecker" and rd != 1:
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
                            inferred_matches = create_pseudobulk_submatrix_vireo(matrix)
                        else:
                            print(f"Error! Do not recognize tool {tool}")
                            continue
                    except FileNotFoundError:
                        print(
                            f"Could not find data for {tool}, pseudobulk dataset {dataset}, ncells {ncells}, read depth {rd}. "
                        )
                        continue

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
