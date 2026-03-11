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
        0,
        index=matrix.index,
        columns=matrix.columns,
    )
    for idx in matrix.index:
        for col in matrix.columns:
            if idx[:-12] == col[:-12]:
                true_matches.loc[idx, col] = 1

    return true_matches


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
        file_path = os.path.join(DATA_PATH + f"2b_bamixchecker/{dataset}/pseudobulk/ncells_{ncells}/BAMixChecker/Total_result.txt")
    else:
        file_path = os.path.join(DATA_PATH + f"2b_bamixchecker/{dataset}/real_data/ncells_null/BAMixChecker/Total_result.txt")

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
            + f"2a_ngscheckmate/{dataset}/real_data/ncells_{ncells}/read_depth_{rd}/output_all.txt",
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
            + f"2a_vireo/{dataset}/real_data/ncells_{ncells}/read_depth_{rd}/similarity_matrix.csv",
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


######################################################
### Functions to calculate accuracy metrics
######################################################
def create_pseudobulk_submatrix_vireo(matrix, hysys_matrix, experiment, n_samples_double=0):
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
    # Filter Vireo matrix to only samples ending in _1 in rows and samples ending in _2 in columns
    # For double samples experiments we want _1 in rows and _2 and _3 in columns.
    if (experiment == "double_samples") and (n_samples_double > 0):
        # For double_samples experiments we have to infer which _3 were included in the experiment, 
        # since Vireo only outputs results for the samples it matched.
        # We do this by looking at the samples included in another tool's results.
        hysys_samples = list(set(hysys_matrix.index))
        hysys_3_samples = [sample for sample in hysys_samples if sample.endswith("_3")]
        matrix_filtered = matrix.loc[
            [idx for idx in matrix.index if idx.endswith("_1")],
            [col for col in matrix.columns if (col.endswith("_2") or col in hysys_3_samples)],
        ]
    else:
        matrix_filtered = matrix.loc[
            [idx for idx in matrix.index if idx.endswith("_1")],
            [col for col in matrix.columns if col.endswith("_2")],
        ]

    # Minimize diagonal of matrix_filtered by moving rows and columns around
    # This is the algorith Vireo uses to match samples
    maximize = False
    if matrix.equals(hysys_matrix):
        maximize = True
    idx0, idx1 = linear_sum_assignment(matrix_filtered.values, maximize)
    matrix_reordered = pd.DataFrame(
        matrix_filtered.iloc[idx0, idx1],
        index=matrix_filtered.index[idx0],
        columns=matrix_filtered.columns[idx1],
    )
    # If there are more columsn than rows, we need to add back the unmatched columns (these will be inferred as non-matches)
    if len(matrix_filtered.columns) > len(matrix_filtered.index):
        # Add back the unmatched columns with NaNs
        unmatched_cols = [col for col in matrix_filtered.columns if col not in matrix_reordered.columns]
        for col in unmatched_cols:
            matrix_reordered[col] = np.nan
    elif len(matrix_filtered.index) > len(matrix_filtered.columns):
        # Add back the unmatched rows with NaNs
        unmatched_rows = [idx for idx in matrix_filtered.index if idx not in matrix_reordered.index]
        for idx in unmatched_rows:
            matrix_reordered.loc[idx] = np.nan
    # The new order of samples in the index and columns reveals the matching of samples between the two pseudobulk sets.
    # Create a sample-matching matrix with 1 on the diagonal and 0 elsewhere, where the order of rows and columns is the same as in matrix_reordered.
    inferred_matches = inferred_matches_pseudobulk_vireo(matrix_reordered)
    # Reorder to the original order of samples
    ordered_columns = sorted(inferred_matches.columns)
    ordered_index = sorted(inferred_matches.index)
    inferred_matches = inferred_matches.loc[ordered_index, ordered_columns]

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


def loop_accuracy_calculations(DATA_PATH, experiment, tools, datasets, ncells_list, read_depths):
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
                            inferred_matches = create_pseudobulk_submatrix_vireo(matrix, hysys_matrix, experiment)
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
                            [idx for idx in inferred_matches.index if idx.endswith("_1")],
                            [col for col in inferred_matches.columns if (col.endswith("_2") or col.endswith("_3"))],
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
    and sort the dataframe by the single cell sample column and whether all rows have data in all columns.

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

    # Find the columns that has "single"
    expected_matches.columns = expected_matches.columns.str.lower()
    if dataset != "hgsoc-new":
        single_cols = [col for col in expected_matches.columns if "single" in col]
    else:
        single_cols = [col for col in expected_matches.columns if "diss" in col]

    # Sort the dataframe first by whether all rows have data in all columns
    # and then by the single cell sample column (assuming there's only one)
    expected_matches["not_all_cols_filled"] = expected_matches.isna().any(axis=1)
    expected_matches = (
        expected_matches.sort_values(by=["not_all_cols_filled", single_cols[0]])
        .drop_duplicates()
        .reset_index(drop=True)
        .drop(columns=["not_all_cols_filled"])
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
        experiment="missing_samples"
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
        print(f"Error! Experiment type {experiment} not recognized. Choose either 'missing_samples' or 'double_samples'.")
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
                            DATA_PATH + f"{folder}{n}/it_{i}/", True, dataset, ncells, rd
                        )
                    )
                elif tool == "HYSYS":
                    inferred_matches = parse_sample_matching_results_hysys(
                        DATA_PATH + f"{folder}{n}/it_{i}/", True, dataset, ncells, rd
                    )
                elif tool == "NGSCheckmate":
                    inferred_matches = (
                        parse_sample_matching_results_ngscheckmate(
                            DATA_PATH + f"{folder}{n}/it_{i}/", True, dataset, ncells, rd
                        )
                    )
                elif tool == "Vireo":
                    matrix = parse_heatmap_matrix_vireo(
                        DATA_PATH + f"{folder}{n}/it_{i}/", True, dataset, ncells, rd
                    )
                    _, hysys_matrix = parse_heatmap_matrix_hysys(
                        DATA_PATH + f"{folder}{n}/it_{i}/", True, dataset, ncells, rd
                    )
                    n_samples_double = n if experiment == "double_samples" else 0
                    inferred_matches = create_pseudobulk_submatrix_vireo(matrix, hysys_matrix, experiment, n_samples_double)
                else:
                    print(
                        f"Error! Do not recognize tool {tool}. Choose one of CrosscheckFingerprints, HYSYS, NGSCheckmate, or Vireo."
                    )
                    continue
                    
                if experiment == "double_samples":
                    inferred_matches_filtered = inferred_matches.loc[
                        [idx for idx in inferred_matches.index if idx.endswith("_1")],
                        [col for col in inferred_matches.columns if (col.endswith("_2") or col.endswith("_3"))],
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

                accuracy_df = pd.concat([accuracy_df,
                    pd.DataFrame({
                        "iteration": [i],
                        "fraction_inconclusive": [fraction_inconclusive],
                        "accuracy": [accuracy],
                        "balanced_accuracy": [balanced_accuracy],
                        "precision": [precision],
                        "recall": [recall],
                        "f1": [f1],
                    })]
                )
    
            final_df = pd.concat([final_df, pd.DataFrame(
                {
                    "tool": [tool],
                    "dataset": [dataset],
                    "ncells": [ncells],
                    "rd": [rd],
                    "n_samples_removed": [n],
                    "av_fraction_inconclusive": accuracy_df["fraction_inconclusive"].mean(),
                    "av_accuracy": accuracy_df["accuracy"].mean(),
                    "av_balanced_accuracy": accuracy_df["balanced_accuracy"].mean(),
                    "av_precision": accuracy_df["precision"].mean(),
                    "av_recall": accuracy_df["recall"].mean(),
                    "av_f1": accuracy_df["f1"].mean(),
                    "sd_fraction_inconclusive": accuracy_df["fraction_inconclusive"].std(),
                    "sd_accuracy": accuracy_df["accuracy"].std(),
                    "sd_balanced_accuracy": accuracy_df["balanced_accuracy"].std(),
                    "sd_precision": accuracy_df["precision"].std(),
                    "sd_recall": accuracy_df["recall"].std(),
                    "sd_f1": accuracy_df["f1"].std(),
                }
            )])
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
    inferred_matches = create_pseudobulk_submatrix_vireo(hysys_matrix, hysys_matrix, "Basic")

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