import os
import re
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

dataset_regex_dict = {
    "hgsoc": r"[a-zA-Z0-9]{32}_",
    "hgsoc-new": r"ds.[a-zA-Z0-9]{32}_",
    "low_grade_glioma": r"^GSM[0-9]{7}_",
}


####################################################
### Functions to convert long dataframes to matrices
####################################################
def long_df_to_matrix_bamixchecker(df, metric="Concordance Rate"):
    """
    Convert long BAMixChecker dataframe to a square matrix format for heatmap visualization.

    input:
        - df: long format dataframe with columns [index (sample), Sample, Concordance Rate, Conclusion]
        - metric: which column to use for the values in the matrix
            (default "Concordance Rate", can also be "Conclusion" for match/no match)

    output:
        - matrix: square dataframe with samples as rows and columns, values are the specified metric
    """
    matrix = df.pivot_table(
        index=0,
        columns="Sample",
        values=metric,
        aggfunc="first",
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
    """
    Convert long CrosscheckFingerprints dataframe to a square matrix format for heatmap visualization.

    input:
        - df: long format dataframe with columns [LEFT_SAMPLE, RIGHT_SAMPLE, LOD_SCORE, RESULT]
        - metric: which column to use for the values in the matrix
            (default "LOD_SCORE", can also be "RESULT" for match/no match)

    output:
        - matrix: square dataframe with samples as rows and columns, values are the specified metric
    """
    matrix = df.pivot_table(
        index="LEFT_SAMPLE",
        columns="RIGHT_SAMPLE",
        values=metric,
        aggfunc="first",
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
    """
    Convert long HYSYS dataframe to a square matrix format for heatmap visualization.

    input:
        - df: long format dataframe with columns [index (sample), Sample, Concordance]

    output:
        - matrix: square dataframe with samples as rows and columns, values are the concordance values
    """
    matrix = df.pivot_table(
        index=0,
        columns="Sample",
        values="Concordance",
        aggfunc="first",
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
    """
    Convert long NGSCheckMate dataframe to a square matrix format for heatmap visualization.

    input:
        - df: long format dataframe with columns [index (sample), Sample, Correlation, Matched]
        - metric: which column to use for the values in the matrix
            (default "Correlation", can also be "Matched" for match/no match)

    output:
        - matrix: square dataframe with samples as rows and columns, values are the specified metric
    """
    matrix = df.pivot_table(
        index=0,
        columns="Sample",
        values=metric,
        aggfunc="first",
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


######################################################
### Function to create true matches matrix for pseudobulk datasets
#######################################################
def true_matches_pseudobulk(matrix):
    """
    Create a true matches (ground truth) matrix for pseudobulk datasets
    based on the pseudobulk names.

    If the samples have identical names apart from the last two characters, then they are a match.

    input:
        - matrix: dataframe with samples as rows and columns

    output:
        - true_matches: dataframe with 1 for true matches and 0 for non-matches
    """
    true_matches = pd.DataFrame(
        np.zeros(matrix.shape),
        index=matrix.index,
        columns=matrix.columns,
    )
    for idx in matrix.index:
        for col in matrix.columns:
            # Get the string before the first underscore in idx
            idx_sample = idx.split("_")[0]
            # Get the string before the first underscore in col
            col_sample = col.split("_")[0]
            if idx_sample == col_sample:
                true_matches.loc[idx, col] = 1

    return true_matches


######################################################
### Functions to help with analysis of Vireo results for pseudobulk datasets
#######################################################
def inferred_matches_pseudobulk_vireo(matrix):
    """
    Create a matrix of inferred matches for Vireo, pseudobulk datasets
    based on the ordering in the heatmap matrix.

    If the number of rows and columns is equal it will simply be an identity matrix.
    If the number of columns is greater than the number of rows, we will add columns of zeroes for the unmatched samples.
    If the number of rows is greater than the number of columns, we will add rows of zeroes for the unmatched samples.

    input:
        - matrix: dataframe with samples as rows and columns

    output:
        - inferred_matches: dataframe with 1 for inferred matches and 0 for non-matches
    """
    min_matrix_size = min(len(matrix.index), len(matrix.columns))
    # Start with a identity matrix
    inferred_matches = np.identity(min_matrix_size)
    inferred_matches = pd.DataFrame(
        inferred_matches,
        index=matrix.index[:min_matrix_size],
        columns=matrix.columns[:min_matrix_size],
    )
    if len(matrix.columns) > len(matrix.index):
        # Add columns with 0s for any samples that don't have a match
        for col in matrix.columns:
            if col not in inferred_matches.columns:
                inferred_matches[col] = 0
    elif len(matrix.index) > len(matrix.columns):
        # Add rows with 0s for any samples that don't have a match
        for idx in matrix.index:
            if idx not in inferred_matches.index:
                inferred_matches.loc[idx] = 0

    return inferred_matches


def create_pseudobulk_submatrix_vireo(
    matrix, hysys_matrix, experiment, n_samples_double=0
):
    """
    Filter Vireo's heatmap matrix to only include comparisons between pseudobulk samples ending in _1
    and samples ending in _2, then apply the same algorithm Vireo uses to infer sample matches
    (minimize diagonal of the matrix) to create a sample-matching matrix with 1 for inferred matches and 0 for non-matches.

    input:
        - matrix: dataframe with samples as rows and columns, values are the similarity scores
        - hysys_matrix: dataframe with samples as rows and columns, values are the similarity scores from HYSYS
        - experiment: the type of experiment (e.g., "double_samples")
        - n_samples_double: the number of double samples in the experiment
    output:
        - inferred_matches: dataframe with 1 for inferred matches and 0 for non-matches
    """
    # Get all unique samples from the HYSYS matrix to filter the Vireo matrix,
    # since Vireo only outputs results for the samples it matched.
    all_samples = list(set(hysys_matrix.index))
    samples_1 = [sample for sample in all_samples if sample.endswith("_1")]
    samples_2 = [sample for sample in all_samples if sample.endswith("_2")]

    # Filter Vireo matrix to only samples ending in _1 in rows and samples ending in _2 in columns
    # For double samples experiments we want _1 in rows and _2 and _3 in columns.
    if (experiment == "double_samples") and (n_samples_double > 0):
        # For double_samples experiments we have to infer which _3 were included in the experiment,
        # since Vireo only outputs results for the samples it matched.
        # We do this by looking at the samples included in another tool's results.
        samples_3 = [sample for sample in all_samples if sample.endswith("_3")]
        matrix_filtered = matrix.loc[
            [idx for idx in matrix.index if idx in samples_1],
            [col for col in matrix.columns if col in samples_2 or col in samples_3],
        ]
        hysys_matrix_filtered = hysys_matrix.loc[
            [idx for idx in hysys_matrix.index if idx in samples_1],
            [
                col
                for col in hysys_matrix.columns
                if col in samples_2 or col in samples_3
            ],
        ]
    elif experiment == "pseudobulk_vs_sc":
        # For pseudobulk vs single-cell experiments we want _1 in rows and _single-cell samples in columns
        matrix_filtered = matrix.loc[
            [idx for idx in matrix.index if idx in samples_1],
            [col for col in matrix.columns if col.endswith("_single-cell")],
        ]
        hysys_matrix_filtered = hysys_matrix.loc[
            [idx for idx in hysys_matrix.index if idx in samples_1],
            [col for col in hysys_matrix.columns if col.endswith("_single-cell")],
        ]
    else:
        matrix_filtered = matrix.loc[
            [idx for idx in matrix.index if idx in samples_1],
            [col for col in matrix.columns if col in samples_2],
        ]
        hysys_matrix_filtered = hysys_matrix.loc[
            [idx for idx in hysys_matrix.index if idx in samples_1],
            [col for col in hysys_matrix.columns if col in samples_2],
        ]

    # Minimize diagonal of matrix_filtered by moving rows and columns around
    # This is the algorith Vireo uses to match samples
    maximize = False
    # For the HYSYS matrix & Vireo algorithm experiment we want to maximize the diagonal instead of minimizing it,
    # since higher concordance values indicate better matches for HYSYS.
    if matrix.equals(hysys_matrix):
        maximize = True
    idx0, idx1 = linear_sum_assignment(matrix_filtered.values, maximize)
    matrix_reordered = pd.DataFrame(
        matrix_filtered.iloc[idx0, idx1],
        index=matrix_filtered.index[idx0],
        columns=matrix_filtered.columns[idx1],
    )
    # If there are more columns than rows (double samples),
    # we need to add back the unmatched columns (these were inferred as non-matches)
    if len(hysys_matrix_filtered.columns) > len(hysys_matrix_filtered.index):
        unmatched_cols = [
            col for col in samples_2 + samples_3 if col not in matrix_reordered.columns
        ]
        for col in unmatched_cols:
            matrix_reordered[col] = 0
    # If there are more rows than columns (missing samples),
    # we need to add back the unmatched rows (these were inferred as non-matches)
    elif len(hysys_matrix_filtered.index) > len(hysys_matrix_filtered.columns):
        unmatched_rows = [idx for idx in samples_1 if idx not in matrix_reordered.index]
        for idx in unmatched_rows:
            matrix_reordered.loc[idx] = 0
    # The new order of samples in the index and columns reveals the matching of samples between the two pseudobulk sets.
    # Create a sample-matching matrix with 1 on the diagonal and 0 elsewhere, where the order of rows and columns is the same as in matrix_reordered.
    inferred_matches = inferred_matches_pseudobulk_vireo(matrix_reordered)
    # Sort the inferred matches
    ordered_columns = sorted(inferred_matches.columns)
    ordered_index = sorted(inferred_matches.index)
    inferred_matches = inferred_matches.loc[ordered_index, ordered_columns]

    return inferred_matches


#####################################################
### Functions to parse tool outputs and create heatmap matrices
#####################################################
def parse_heatmap_matrix_bamixchecker(DATA_PATH, pseudobulk, dataset, ncells):
    """
    Read BAMixChecker output and create a matrix for heatmap visualization.

    input:
        - DATA_PATH: base path to the data directory
        - pseudobulk: boolean indicating if the data is pseudobulk
        - dataset: the dataset name
        - ncells: the number of cells used to create the pseudobulk (or "null" for real data)

    output:
        - df: long format dataframe with columns [index (sample), Sample, Concordance Rate, Conclusion]
        - matrix: square dataframe with samples as rows and columns, values are the concordance rates
    """
    if pseudobulk:
        file_path = os.path.join(
            DATA_PATH,
            f"2b_bamixchecker/{dataset}/pseudobulk/ncells_{ncells}/BAMixChecker/Total_result.txt",
        )
    else:
        file_path = os.path.join(
            DATA_PATH,
            f"2b_bamixchecker/{dataset}/real_data/ncells_null/BAMixChecker/Total_result.txt",
        )

    df = pd.read_csv(
        file_path,
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
    """
    Read CrosscheckFingerprints output and create a matrix for heatmap visualization.

    input:
        - DATA_PATH: base path to the data directory
        - pseudobulk: boolean indicating if the data is pseudobulk
        - dataset: the dataset name
        - ncells: the number of cells used to create the pseudobulk (or "null" for real data)
        - rd: the read depth filter cut-off used for the analysis

    output:
        - df: long format dataframe with columns [LEFT_SAMPLE, RIGHT_SAMPLE, LOD_SCORE, RESULT]
        - matrix: square dataframe with samples as rows and columns, values are the LOD_SCOREs
    """
    if pseudobulk:
        file_path = os.path.join(
            DATA_PATH
            + f"2a_fingerprints/{dataset}/pseudobulk/ncells_{ncells}/read_depth_{rd}/crosscheck_metrics.txt"
        )
    else:
        file_path = os.path.join(
            DATA_PATH
            + f"2a_fingerprints/{dataset}/real_data/ncells_null/read_depth_{rd}/crosscheck_metrics.txt"
        )
    df = pd.read_csv(
        file_path,
        sep="\t",
        skiprows=6,
        header=0,
    )

    df = df[["LEFT_SAMPLE", "RIGHT_SAMPLE", "LOD_SCORE", "RESULT"]]

    if pseudobulk or dataset == "hgsoc-new":
        regex_exp = dataset_regex_dict.get(dataset, "")
    else:
        regex_exp = ""

    rename_dict = {
        x: re.sub(regex_exp, "", x.replace(".bam", ""))
        for x in list(set(df["LEFT_SAMPLE"]) | set(df["RIGHT_SAMPLE"]))
    }
    df["LEFT_SAMPLE"] = df["LEFT_SAMPLE"].map(rename_dict)
    df["RIGHT_SAMPLE"] = df["RIGHT_SAMPLE"].map(rename_dict)

    # Create matrix for heatmap
    matrix = long_df_to_matrix_crosscheckfingerprints(df, metric="LOD_SCORE")

    return df, matrix


def parse_heatmap_matrix_hysys(DATA_PATH, pseudobulk, dataset, ncells, rd):
    """
    Read HYSYS output and create a matrix for heatmap visualization.

    input:
        - DATA_PATH: base path to the data directory
        - pseudobulk: boolean indicating if the data is pseudobulk
        - dataset: the dataset name
        - ncells: the number of cells used to create the pseudobulk (or "null" for real data)
        - rd: the read depth filter cut-off used for the analysis

    output:
        - df: long format dataframe with columns [Sample, Concordance, Size]
        - matrix: square dataframe with samples as rows and columns, values are the Concordance values
    """
    if pseudobulk:
        file_path = os.path.join(
            DATA_PATH
            + f"2a_hysys/{dataset}/pseudobulk/ncells_{ncells}/read_depth_{rd}/concordance_output.txt"
        )
    else:
        file_path = os.path.join(
            DATA_PATH
            + f"2a_hysys/{dataset}/real_data/ncells_null/read_depth_{rd}/concordance_output.txt"
        )

    df = pd.read_csv(
        file_path,
        sep="\t",
        index_col=0,
        header=None,
    )

    if pseudobulk or dataset == "hgsoc-new":
        regex_exp = dataset_regex_dict.get(dataset, "")
    else:
        regex_exp = ""

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
    """
    Read NGSCheckmate output and create a matrix for heatmap visualization.

    input:
        - DATA_PATH: base path to the data directory
        - pseudobulk: boolean indicating if the data is pseudobulk
        - dataset: the dataset name
        - ncells: the number of cells used to create the pseudobulk (or "null" for real data)
        - rd: the read depth filter cut-off used for the analysis

    output:
        - df: long format dataframe with columns [Matched, Sample, Binary, Correlation]
        - matrix: square dataframe with samples as rows and columns, values are the Correlation values
    """
    if pseudobulk:
        file_path = os.path.join(
            DATA_PATH
            + f"2a_ngscheckmate/{dataset}/pseudobulk/ncells_{ncells}/read_depth_{rd}/output_all.txt",
        )
    else:
        file_path = os.path.join(
            DATA_PATH
            + f"2a_ngscheckmate/{dataset}/real_data/ncells_null/read_depth_{rd}/output_all.txt",
        )

    df = pd.read_csv(
        file_path,
        sep="\t",
        index_col=0,
        header=None,
    )

    if pseudobulk or dataset == "hgsoc-new":
        regex_exp = dataset_regex_dict.get(dataset, "")
    else:
        regex_exp = ""

    # Clean up sample names
    rename_dict = {
        x: re.sub(
            regex_exp,
            "",
            x.replace(".vcf", "").replace("_individual_variants", ""),
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
    """
    Read Vireo output and create a matrix for heatmap visualization.

    input:
        - DATA_PATH: base path to the data directory
        - pseudobulk: boolean indicating if the data is pseudobulk
        - dataset: the dataset name
        - ncells: the number of cells used to create the pseudobulk (or "null" for real data)
        - rd: the read depth filter cut-off used for the analysis

    output:
        - df: long format dataframe with columns [sample_id_0, sample_id_1, match]
        - matrix: dataframe with samples from one modality as rows and samples from another modality as columns,
            values are the similarity scores used by Vireo to match samples (lower values indicate more similar samples)
    """
    if pseudobulk:
        file_path = os.path.join(
            DATA_PATH
            + f"2a_vireo/{dataset}/pseudobulk/ncells_{ncells}/read_depth_{rd}/similarity_matrix.csv",
        )
    else:
        file_path = os.path.join(
            DATA_PATH
            + f"2a_vireo/{dataset}/real_data/ncells_null/read_depth_{rd}/similarity_matrix.csv",
        )
    matrix = pd.read_csv(
        file_path,
        index_col=0,
    )

    if pseudobulk or dataset == "hgsoc-new":
        regex_exp = dataset_regex_dict.get(dataset, "")
    else:
        regex_exp = ""

    matrix.columns = [
        re.sub(regex_exp, "", col.replace(".bam", "")) for col in matrix.columns
    ]
    matrix.index = pd.Index(
        [re.sub(regex_exp, "", idx.replace(".bam", "")) for idx in matrix.index]
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
    """
    Parse BAMixChecker Total_result.txt file to categorize sample relationships.

    input:
        - DATA_PATH: base path to the data directory
        - pseudobulk: boolean indicating if the data is pseudobulk
        - dataset: the dataset name
        - ncells: the number of cells used to create the pseudobulk (or "null" for real data)

    output:
        - df: square dataframe with binary values indicating sample matches (1 for match, 0 for no match)
    """
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
    """
    Parse CrosscheckFingerprints results to categorize sample relationships.

    input:
        - DATA_PATH: base path to the data directory
        - pseudobulk: boolean indicating if the data is pseudobulk
        - dataset: the dataset name
        - ncells: the number of cells used to create the pseudobulk (or "null" for real data)

    output:
        - df: square dataframe with binary values indicating sample matches (1 for match, 0 for no match)
    """
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

    input:
        - DATA_PATH: base path to the data directory
        - pseudobulk: boolean indicating if the data is pseudobulk
        - dataset: the dataset name
        - ncells: the number of cells used to create the pseudobulk (or "null" for real data)

    output:
        - df: square dataframe with binary values indicating sample matches (1 for match, 0 for no match)
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
        if pair in matching_pairs or pair[::-1] in matching_pairs:
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
    """
    Parse NGSCheckMate results to categorize sample relationships.

    input:
        - DATA_PATH: base path to the data directory
        - pseudobulk: boolean indicating if the data is pseudobulk
        - dataset: the dataset name
        - ncells: the number of cells used to create the pseudobulk (or "null" for real data)

    output:
        - df: square dataframe with binary values indicating sample matches (1 for match, 0 for no match)
    """

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
    """
    Parse Vireo's matched_samples.csv to categorize sample relationships.

    input:
        - DATA_PATH: base path to the data directory
        - pseudobulk: boolean indicating if the data is pseudobulk
        - dataset: the dataset name
        - ncells: the number of cells used to create the pseudobulk (or "null" for real data)

    output:
        - df: square dataframe with binary values indicating sample matches (1 for match, 0 for no match)
    """
    if pseudobulk:
        file_path = os.path.join(
            DATA_PATH
            + f"2a_vireo/{dataset}/pseudobulk/ncells_{ncells}/read_depth_{rd}/matched_samples.csv",
        )
    else:
        file_path = os.path.join(
            DATA_PATH
            + f"2a_vireo/{dataset}/real_data/ncells_{ncells}/read_depth_{rd}/matched_samples.csv",
        )
    sample_matches = pd.read_csv(
        file_path,
        index_col=0,
    ).reset_index()

    if pseudobulk:
        regex_exp = dataset_regex_dict.get(dataset, "")
    else:
        regex_exp = ""
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
