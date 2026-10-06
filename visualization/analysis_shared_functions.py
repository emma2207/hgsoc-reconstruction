import os
import re
import glob
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

dataset_regex_dict = {
    "hgsoc": r"[a-zA-Z0-9]{32}_",
    "hgsoc-new": r"ds.[a-zA-Z0-9]{32}_",
    # "low_grade_glioma": r"^GSM[0-9]{7}_",
}


def _analysis_base_path(dataset, pseudobulk, ncells, mod1, mod2):
    """Return the default dataset-relative base folder used by analysis outputs."""
    if pseudobulk:
        return os.path.join(dataset, "pseudobulk", f"ncells_{ncells}")
    return os.path.join(dataset, f"{mod1}_vs_{mod2}", "ncells_null")


def _modality_path_aliases(modality):
    aliases = {
        "bulk_diss_polyA": "bulk_dissociated_polyA",
        "bulk_dissociated_polyA": "bulk_diss_polyA",
    }
    return [modality, aliases[modality]] if modality in aliases else [modality]


def _candidate_analysis_base_paths(DATA_PATH, tool, dataset, pseudobulk, ncells, mod1, mod2):
    """Return candidate base folders in priority order for analysis outputs."""
    if pseudobulk:
        return [
            os.path.join(
                DATA_PATH,
                tool,
                _analysis_base_path(dataset, pseudobulk, ncells, mod1, mod2),
            )
        ]

    candidate_rel_paths = []
    for left in _modality_path_aliases(mod1):
        for right in _modality_path_aliases(mod2):
            candidate_rel_paths.extend(
                [
                    os.path.join(dataset, f"{left}_vs_{right}", "ncells_null"),
                    os.path.join(dataset, f"{right}_vs_{left}", "ncells_null"),
                ]
            )
    candidate_paths = [
        os.path.join(DATA_PATH, tool, rel_path) for rel_path in candidate_rel_paths
    ]

    # Keep order while removing duplicates.
    deduped_paths = []
    for path in candidate_paths:
        if path not in deduped_paths:
            deduped_paths.append(path)

    return deduped_paths


def _resolve_analysis_base_path(DATA_PATH, tool, dataset, pseudobulk, ncells, mod1, mod2):
    """Return the first existing tool output base path for the requested modalities."""
    candidate_paths = _candidate_analysis_base_paths(
        DATA_PATH, tool, dataset, pseudobulk, ncells, mod1, mod2
    )
    return next((path for path in candidate_paths if os.path.isdir(path)), candidate_paths[0])


def _normalize_read_depth_tag(rd):
    """Normalize read-depth value to a tag without the read_depth_ prefix."""
    rd_tag = str(rd)
    if rd_tag.startswith("read_depth_"):
        rd_tag = rd_tag.replace("read_depth_", "", 1)
    return rd_tag


def _modality_in_label(modality, label):
    """Return True when a sample label appears to belong to the requested modality."""
    modality_key = str(modality).lower().replace("-", "_").replace(" ", "_")
    label_key = str(label).lower().replace("-", "_").replace(" ", "_")

    if modality_key in label_key:
        return True

    # Treat single-cell/nucleus aliases as compatible, including pooled_single_cell.
    modality_tokens = modality_key.split("_")
    if "single" in modality_tokens and "single" in label_key:
        return True

    return False


def _resolve_read_depth_path(
    DATA_PATH, tool, dataset, pseudobulk, ncells, rd, mod1, mod2
):
    """Resolve a read-depth directory for both legacy scalar and new modality-specific tags."""
    rd_tag = _normalize_read_depth_tag(rd)
    candidate_base_paths = _candidate_analysis_base_paths(
        DATA_PATH,
        tool,
        dataset,
        pseudobulk,
        ncells,
        mod1,
        mod2,
    )

    for base_path in candidate_base_paths:
        if not os.path.exists(base_path):
            continue

        prefixed_matches = sorted(
            glob.glob(os.path.join(base_path, f"read_depth_{rd_tag}*"))
        )
        if len(prefixed_matches) == 1:
            return prefixed_matches[0]

        # A scalar read depth denotes that depth for both requested modalities.
        # Prefer the exact modality-qualified folder before trying ambiguous token matches.
        if rd_tag.isdigit():
            exact_depth_paths = [
                os.path.join(
                    base_path,
                    f"read_depth_{left}_{rd_tag}_{right}_{rd_tag}",
                )
                for left in _modality_path_aliases(mod1)
                for right in _modality_path_aliases(mod2)
            ]
            exact_depth_matches = [
                path for path in exact_depth_paths if os.path.isdir(path)
            ]
            if len(exact_depth_matches) == 1:
                return exact_depth_matches[0]

        all_matches = sorted(glob.glob(os.path.join(base_path, "read_depth_*")))

        # Support modality-specific tags when rd is scalar, e.g.
        # read_depth_bulk_diss_polyA_20_single-cell_20 for rd=20.
        token_matches = []
        for candidate in all_matches:
            candidate_tag = os.path.basename(candidate).replace("read_depth_", "", 1)
            tokens = re.split(r"[_-]", candidate_tag)
            if rd_tag in tokens:
                token_matches.append(candidate)
        if len(token_matches) == 1:
            return token_matches[0]

        # If modality labels differ from folder naming (e.g. single-cell vs single-nucleus),
        # fall back to matching only by numeric depth tokens.
        rd_numbers = re.findall(r"\d+", rd_tag)
        if rd_numbers:
            numeric_matches = []
            for candidate in all_matches:
                candidate_tag = os.path.basename(candidate).replace(
                    "read_depth_", "", 1
                )
                candidate_numbers = re.findall(r"\d+", candidate_tag)
                if candidate_numbers == rd_numbers:
                    numeric_matches.append(candidate)
            if len(numeric_matches) == 1:
                return numeric_matches[0]

        if len(all_matches) == 1:
            return all_matches[0]

    candidate_base_str = ", ".join(candidate_base_paths)
    raise FileNotFoundError(
        f"Could not resolve read depth path for rd={rd}. Checked base paths: {candidate_base_str}."
    )


def _split_sample_label(label):
    """Split sample label into prefix/suffix at first underscore, safely."""
    label = str(label)
    parts = label.split("_", 1)
    if len(parts) == 2:
        return parts[0], parts[1]
    return label, ""


def _sorted_sample_labels(labels):
    """Sort labels by suffix then prefix; tolerate labels without underscores."""
    return sorted(
        labels,
        key=lambda x: (
            _split_sample_label(x)[1],
            _split_sample_label(x)[0],
            str(x),
        ),
    )


def matches_matrix_to_pair_df(matrix):
    """
    Convert a 0/1 sample-by-sample match matrix into a long dataframe of matched sample pairs.

    input:
        - matrix: dataframe where rows/columns are sample IDs and values are 0/1

    output:
        - matches_df: dataframe with one row per inferred match
            If matrix index/columns share names, those are used as output column names.
            Otherwise defaults to ["sample_id_0", "sample_id_1"].
    """
    matrix_df = matrix.copy()

    matrix_df = matrix_df.apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)

    row_name = matrix_df.index.name or "sample_id_0"
    col_name = matrix_df.columns.name or "sample_id_1"

    matches_df = matrix_df.eq(1).stack().loc[lambda s: s].index.to_frame(index=False)
    matches_df.columns = [row_name, col_name]

    if matrix_df.index.equals(matrix_df.columns):
        matches_df = matches_df[
            matches_df[row_name] < matches_df[col_name]
        ].reset_index(drop=True)

    return matches_df


####################################################
### Functions to convert long dataframes to matrices
####################################################
def _long_df_to_matrix(
    df,
    row_column,
    column_column,
    value_column,
    *,
    symmetric=True,
    diagonal_value=None,
    value_scale: float = 1.0,
    lexical_sort=False,
):
    """Build a labeled square matrix from pairwise long-form data."""
    row_labels = df.index if row_column is None else df[row_column]
    column_labels = df[column_column]
    samples = sorted(set(row_labels) | set(column_labels))

    long_data = pd.DataFrame(
        {
            "_row": row_labels.to_numpy(),
            "_column": column_labels.to_numpy(),
            "_value": df[value_column].to_numpy(),
        }
    )
    matrix = long_data.pivot_table(
        index="_row",
        columns="_column",
        values="_value",
        aggfunc="first",
    ).reindex(index=samples, columns=samples)
    if value_scale != 1:
        matrix = matrix * value_scale
    if diagonal_value is not None:
        values = matrix.to_numpy(copy=True)
        np.fill_diagonal(values, diagonal_value)
        matrix = pd.DataFrame(values, index=matrix.index, columns=matrix.columns)
    if symmetric:
        matrix = matrix.combine_first(matrix.T)

    ordered_labels = (
        sorted(matrix.index)
        if lexical_sort
        else _sorted_sample_labels(matrix.index)
    )
    matrix = matrix.reindex(index=ordered_labels, columns=ordered_labels)
    matrix.index.name = df.index.name if row_column is None else row_column
    matrix.columns.name = column_column
    return matrix


def long_df_to_matrix_bamixchecker(df, metric="Concordance Rate"):
    return _long_df_to_matrix(
        df, None, "Sample", metric, lexical_sort=True
    )


def long_df_to_matrix_conpair(df, metric="concordance"):
    return _long_df_to_matrix(df, "sample1", "sample2", metric)


def long_df_to_matrix_crosscheckfingerprints(df, metric="LOD_SCORE"):
    return _long_df_to_matrix(df, "LEFT_SAMPLE", "RIGHT_SAMPLE", metric)


def long_df_to_matrix_hysys(df):
    return _long_df_to_matrix(
        df, None, "Sample", "Concordance"
    )


def long_df_to_matrix_ngscheckmate(df, metric="Correlation"):
    value_scale = 0.01 if metric == "Correlation" else 1
    return _long_df_to_matrix(
        df, None, "Sample", metric, value_scale=value_scale
    )


def long_df_to_matrix_ntsm(df):
    return _long_df_to_matrix(
        df, "sample1", "sample2", "relate", symmetric=False
    )


def long_df_to_matrix_omicsprint(df):
    return _long_df_to_matrix(df, "sample1", "sample2", "mean")


def long_df_to_matrix_peddy(df, metric="rel_difference"):
    diagonal_values = {"rel_difference": 0, "match": 1}
    if metric not in diagonal_values:
        raise ValueError(f"Unsupported Peddy matrix metric: {metric}")
    return _long_df_to_matrix(
        df,
        "sample_a",
        "sample_b",
        metric,
        diagonal_value=diagonal_values[metric],
    )


def long_df_to_matrix_somalier(df, metric="relatedness"):
    return _long_df_to_matrix(
        df,
        "sample_a",
        "sample_b",
        metric,
        diagonal_value=0,
    )


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
def parse_heatmap_matrix_bamixchecker(
    DATA_PATH,
    pseudobulk,
    dataset,
    ncells,
    mod1="bulk_chunk_ribo",
    mod2="bulk_diss_polyA",
):
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
        base_path = _resolve_analysis_base_path(
            DATA_PATH,
            "2b_bamixchecker",
            dataset,
            pseudobulk,
            ncells,
            mod1,
            mod2,
        )
        file_path = os.path.join(base_path, "BAMixChecker", "Total_result.txt")

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


def parse_heatmap_matrix_conpair(
        DATA_PATH,
        pseudobulk,
        dataset,
        ncells,
        mod1="bulk_chunk_ribo",
        mod2="bulk_diss_polyA",
):
    """
    Read Conpair output and create a matrix for heatmap visualization.

    input:
        - DATA_PATH: base path to the data directory
        - pseudobulk: boolean indicating if the data is pseudobulk
        - dataset: the dataset name
        - ncells: the number of cells used to create the pseudobulk (or "null" for real data)
        - mod1: the first modality
        - mod2: the second modality

    output:
        - long_df: long format dataframe with columns [sample1, sample2, concordance]
        - matrix: square dataframe with samples as rows and columns, values are the concordance rates
    """
    
    if pseudobulk:
        file_path = os.path.join(
            DATA_PATH,
            f"2b_conpair/{dataset}/pseudobulk/ncells_{ncells}",
        )
    else:
        file_path = _resolve_analysis_base_path(
            DATA_PATH,
            "2b_conpair",
            dataset,
            pseudobulk,
            ncells,
            mod1,
            mod2,
        )
    long_df = pd.DataFrame()

    # Loop over all files in the directory
    for filename in os.listdir(file_path): 
        if filename.endswith(".txt"):
            # Extract sample names from the filename
            # Pair-index prefix is "pb_<i>_<j>_" for pseudobulk pairs and "rd_<i>_<j>_"
            # for real-data pairs (an arbitrary pairing index, not a read-depth value).
            sample1 = filename.split("_vs_")[0]
            sample1 = re.sub(r"(pb|rd)_\d{1,2}_\d{1,2}_A_", "", sample1).replace("_filtered", "")
            if dataset == "hgsoc-new":
                sample1 = re.sub(r"_ds\.[a-z0-9]{32}", "", sample1)
            sample2 = filename.split("_vs_")[1]
            sample2 = re.sub(r"(pb|rd)_\d{1,2}_\d{1,2}_B_", "", sample2).replace("_filtered_concordance.txt", "")
            if dataset == "hgsoc-new":
                sample2 = re.sub(r"_ds\.[a-z0-9]{32}", "", sample2)

            # Tag real-data sample names with their modality so downstream modality-based
            # filtering/ordering can distinguish mod1 (A) from mod2 (B) samples.
            # Pseudobulk sample names keep their "_1"/"_2" suffix untouched.
            if not pseudobulk:
                def collapse_modality_suffix(sample_name):
                    for modality in (mod1, mod2):
                        variants = list(dict.fromkeys((
                            modality,
                            modality.replace("-", "_"),
                            modality.replace("_", "-"),
                        )))
                        for first in variants:
                            for second in variants:
                                repeated = f"_{first}_{second}"
                                while repeated in sample_name:
                                    sample_name = sample_name.replace(
                                        repeated, f"_{first}"
                                    )
                    return sample_name

                sample1 = collapse_modality_suffix(sample1)
                sample2 = collapse_modality_suffix(sample2)
                modality_suffixes = {
                    f"_{variant}"
                    for modality in (mod1, mod2)
                    for variant in (
                        modality,
                        modality.replace("-", "_"),
                        modality.replace("_", "-"),
                    )
                }
                if not sample1.endswith(tuple(modality_suffixes)):
                    sample1 = f"{sample1}_{mod1}"
                if not sample2.endswith(tuple(modality_suffixes)):
                    sample2 = f"{sample2}_{mod2}"

            # Read the concordance value from the file
            INPUT_FILE = os.path.join(file_path, filename)
            df = pd.read_csv(INPUT_FILE, sep="\t", header=None)
            concordance_string = str(df.iloc[0, 0])
            concordance_string = concordance_string.replace("Concordance: ", "").replace("%", "")
            concordance_row = pd.DataFrame({
                "sample1": [sample1],
                "sample2": [sample2],
                "concordance": [float(concordance_string)],
            })
            long_df = pd.concat([long_df, concordance_row], ignore_index=True)

    matrix = long_df_to_matrix_conpair(long_df)

    return long_df, matrix


def parse_heatmap_matrix_crosscheckfingerprints(
    DATA_PATH,
    pseudobulk,
    dataset,
    ncells,
    rd,
    mod1="bulk_chunk_ribo",
    mod2="bulk_diss_polyA",
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
        rd_path = f"{DATA_PATH}/2a_fingerprints/{dataset}/pseudobulk/ncells_{ncells}/read_depth_{rd}"
    else:
        rd_path = _resolve_read_depth_path(
            DATA_PATH, "2a_fingerprints", dataset, pseudobulk, ncells,
            f"{mod1}_0_{mod2}_0", mod1, mod2,
        )
    file_path = os.path.join(rd_path, "crosscheck_metrics.txt")
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
    df["LOD_SCORE"] = df["LOD_SCORE"].astype(float)
    if dataset == "hgsoc":
        df["LEFT_SAMPLE"] = df["LEFT_SAMPLE"].str.replace(f"_{mod1}", "", 1).str.replace(f"_{mod2}", "", 1)
        df["RIGHT_SAMPLE"] = df["RIGHT_SAMPLE"].str.replace(f"_{mod1}", "", 1).str.replace(f"_{mod2}", "", 1)

    # Create matrix for heatmap
    matrix = long_df_to_matrix_crosscheckfingerprints(df, metric="LOD_SCORE")

    return df, matrix


def parse_heatmap_matrix_hysys(
    DATA_PATH,
    pseudobulk,
    dataset,
    ncells,
    rd,
    mod1="bulk_chunk_ribo",
    mod2="bulk_diss_polyA",
):
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
        rd_path = f"{DATA_PATH}/2a_hysys/{dataset}/pseudobulk/ncells_{ncells}/read_depth_{rd}"
    else:
        rd_path = _resolve_read_depth_path(
            DATA_PATH, "2a_hysys", dataset, pseudobulk, ncells,
            f"{mod1}_0_{mod2}_0", mod1, mod2,
        )
    file_path = os.path.join(rd_path, "concordance_output.txt")

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
            os.path.basename(x).replace("_individual_variants.snps", ""),
        )
        for x in list(set(df.index) | set(df[1]))
    }
    df.index = df.index.map(rename_dict)
    df[1] = df[1].map(rename_dict)
    df.columns = ["Sample", "Concordance", "Size"]

    matrix = long_df_to_matrix_hysys(df)

    return df, matrix


def parse_heatmap_matrix_ngscheckmate(
    DATA_PATH,
    pseudobulk,
    dataset,
    ncells,
    rd,
    mod1="bulk_chunk_ribo",
    mod2="bulk_diss_polyA",
):
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
        rd_path = f"{DATA_PATH}/2a_ngscheckmate/{dataset}/pseudobulk/ncells_{ncells}/read_depth_{rd}"
    else:
        rd_path = _resolve_read_depth_path(
            DATA_PATH, "2a_ngscheckmate", dataset, pseudobulk, ncells,
            f"{mod1}_0_{mod2}_0", mod1, mod2,
        )
    file_path = os.path.join(rd_path, "output_all.txt")
    correlation_matrix_path = os.path.join(
        rd_path, "output_output_corr_matrix.txt"
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

    # Prefer the full-precision correlation matrix over rounded values in output_all.txt.
    matrix = long_df_to_matrix_ngscheckmate(df)
    if os.path.exists(correlation_matrix_path):
        matrix = pd.read_csv(correlation_matrix_path, sep="\t", index_col=0)
        clean_sample_name = lambda sample: re.sub(
            regex_exp,
            "",
            str(sample).replace(".vcf", "").replace("_individual_variants", ""),
        )
        matrix.index = matrix.index.map(clean_sample_name)
        matrix.columns = matrix.columns.map(clean_sample_name)
        ordered_labels = _sorted_sample_labels(matrix.index)
        matrix = matrix.reindex(index=ordered_labels, columns=ordered_labels)

    return df, matrix


def parse_heatmap_matrix_ntsm(
    DATA_PATH,
    pseudobulk,
    dataset,
    ncells,
    mod1="bulk_chunk_ribo",
    mod2="bulk_diss_polyA",
):
    """
    Read NTSM output and create a matrix for heatmap visualization.

    input:
        - DATA_PATH: base path to the data directory
        - pseudobulk: boolean indicating if the data is pseudobulk
        - dataset: the dataset name
        - ncells: the number of cells used to create the pseudobulk (or "null" for real data)
        - mod1: the first modality
        - mod2: the second modality

    output:
        - df: long format dataframe with columns [sample_a, sample_b, score]
        - matrix: square dataframe with samples as rows and columns, values are the loglikelihood scores between samples
    """

    if pseudobulk:
        file_path = os.path.join(
            DATA_PATH,
            f"ntsm_eval/ntsm_pairwise_pseudobulk_{dataset}_ncells_{ncells}.tsv",
        )
    else:
        modality_orders = [(mod1, mod2)]
        if dataset in {"low_grade_glioma", "hgsoc", "hgsoc-new"}:
            modality_orders.insert(0, (mod2, mod1))
        candidate_paths = [
            os.path.join(
                DATA_PATH,
                f"ntsm_eval/ntsm_pairwise_{left}_{right}_{dataset}.tsv",
            )
            for left, right in modality_orders
        ]
        file_path = next(
            (candidate for candidate in candidate_paths if os.path.exists(candidate)),
            candidate_paths[0],
        )
    ntsm_df = pd.read_csv(
        file_path,
        sep="\t",
        header=0,
    )
    ntsm_df["sample1"] = ntsm_df["sample1"].str.replace("counts_", "").str.replace(".txt", "")
    ntsm_df["sample2"] = ntsm_df["sample2"].str.replace("counts_", "").str.replace(".txt", "")
    ntsm_df = ntsm_df[["sample1", "sample2", "relate"]]
    # Symmetrize the dataframe by adding the reverse pairs
    ntsm_df = pd.concat([ntsm_df, ntsm_df.rename(columns={"sample1": "sample2", "sample2": "sample1"})], ignore_index=True)

    def fix_sample_id(sample_id):
        parts = sample_id.split("_")

        if parts[-1].startswith("HGSOC-"):
            return f"{parts[-1]}_{'_'.join(parts[:-1])}"

        return sample_id


    ntsm_df["sample1"] = ntsm_df["sample1"].apply(fix_sample_id)
    ntsm_df["sample2"] = ntsm_df["sample2"].apply(fix_sample_id)

    # ntsm only puts out matching pairs, not non-matching pairs
    # so we need to read the missing sample pairs from a different tool
    crosscheck_df, _ = parse_heatmap_matrix_crosscheckfingerprints(
        DATA_PATH,
        pseudobulk,
        dataset,
        ncells,
        rd=0,
        mod1=mod1,
        mod2=mod2,
    )

    sample_label_by_id = {}
    crosscheck_sample_labels = set(crosscheck_df["LEFT_SAMPLE"]) | set(
        crosscheck_df["RIGHT_SAMPLE"]
    )
    for sample_label in crosscheck_sample_labels:
        for modality in (mod1, mod2):
            suffix = f"_{modality}"
            if sample_label.endswith(suffix):
                sample_label_by_id[sample_label[: -len(suffix)]] = sample_label

    def canonical_sample_label(sample_id):
        sample_id = str(sample_id)
        for modality in (mod1, mod2):
            prefix = f"{modality}_"
            if sample_id.startswith(prefix):
                sample_id = sample_id.removeprefix(prefix)
                break
        return sample_label_by_id.get(sample_id, sample_id)

    ntsm_df["sample1"] = ntsm_df["sample1"].map(canonical_sample_label)
    ntsm_df["sample2"] = ntsm_df["sample2"].map(canonical_sample_label)

    crosscheck_pairs = set()
    for _, row in crosscheck_df.iterrows():
        left, right = row["LEFT_SAMPLE"], row["RIGHT_SAMPLE"]
        crosscheck_pairs.add(tuple(sorted([left, right])))
    # Add missing sample pairs from crosscheck_pairs to the ntsm_df
    existing_pairs = set(
        tuple(sorted([row["sample1"], row["sample2"]])) for _, row in ntsm_df.iterrows()
    )
    for pair in crosscheck_pairs:
        if pair not in existing_pairs:
            ntsm_df = pd.concat(
                [
                    ntsm_df,
                    pd.DataFrame(
                        [
                            {"sample1": pair[0], "sample2": pair[1], "relate": np.nan},
                            {"sample1": pair[1], "sample2": pair[0], "relate": np.nan},
                        ]
                    ),
                ],
                ignore_index=True,
            )
    matrix = long_df_to_matrix_ntsm(ntsm_df)

    return ntsm_df, matrix.astype(float)


def parse_heatmap_matrix_omicsprint(
    DATA_PATH,
    pseudobulk,
    dataset,
    ncells,
    rd,
    mod1="bulk_chunk_ribo",
    mod2="bulk_diss_polyA",
):
    """
    Read OmicsPrint output and create a matrix for heatmap visualization.

    input:
        - DATA_PATH: base path to the data directory
        - pseudobulk: boolean indicating if the data is pseudobulk
        - dataset: the dataset name
        - ncells: the number of cells used to create the pseudobulk (or "null" for real data)
        - rd: the read depth filter cut-off used for the analysis

    output:
        - df: long format dataframe with columns [sample1, sample2, mean, relation]
        - matrix: square dataframe with samples as rows and columns, values are the mean ibs (identity by state) values
    """
    if pseudobulk:
        rd_path = os.path.join(
            DATA_PATH,
            f"2a_omicsprint/{dataset}/pseudobulk/ncells_{ncells}/read_depth_{rd}",
        )
    else:
        rd_path = _resolve_read_depth_path(
            DATA_PATH, "2a_omicsprint", dataset, pseudobulk, ncells,
            f"{mod1}_0_{mod2}_0", mod1, mod2,
        )
    nested_file_path = os.path.join(
        rd_path, "omicsprint", "omicsprint_allele_sharing.tsv"
    )
    flat_file_path = os.path.join(rd_path, "omicsprint_allele_sharing.tsv")
    file_path = nested_file_path if os.path.exists(nested_file_path) else flat_file_path
    df = pd.read_csv(
        file_path,
        sep="\t",
        header=0,
    )

    df = df[["colnames.x", "colnames.y", "mean"]]
    df.columns = ["sample1", "sample2", "mean"]

    if pseudobulk or dataset == "hgsoc-new":
        regex_exp = dataset_regex_dict.get(dataset, "")
    else:
        regex_exp = ""

    rename_dict = {
        x: re.sub(regex_exp, "", x.replace(".bam", ""))
        for x in list(set(df["sample1"]) | set(df["sample2"]))
    }
    df["sample1"] = df["sample1"].map(rename_dict)
    df["sample2"] = df["sample2"].map(rename_dict)

    if dataset == "hgsoc":
        df["sample1"] = df["sample1"].str.replace(f"_{mod1}", "", 1).str.replace(f"_{mod2}", "", 1)
        df["sample2"] = df["sample2"].str.replace(f"_{mod1}", "", 1).str.replace(f"_{mod2}", "", 1)

    # Get all samples from CrosscheckFingerprints output to ensure that the OmicsPrint matrix has the same samples as the other tools.
    crosscheck_df, _ = parse_heatmap_matrix_crosscheckfingerprints(
        DATA_PATH,
        pseudobulk,
        dataset,
        ncells,
        rd,
        mod1=mod1,
        mod2=mod2,
    )
    crosscheck_pairs = set(tuple(sorted([row["LEFT_SAMPLE"], row["RIGHT_SAMPLE"]])) for _, row in crosscheck_df.iterrows())

    # Add missing sample pairs to the omicsprint df
    for pair in crosscheck_pairs:
        if pair not in set(tuple(sorted([row["sample1"], row["sample2"]])) for _, row in df.iterrows()):
            df = pd.concat([df, pd.DataFrame([{"sample1": pair[0], "sample2": pair[1], "mean": np.nan}])], ignore_index=True)
            df = pd.concat([df, pd.DataFrame([{"sample1": pair[1], "sample2": pair[0], "mean": np.nan}])], ignore_index=True)

    # Create matrix for heatmap
    matrix = long_df_to_matrix_omicsprint(df)

    return df, matrix


def parse_heatmap_matrix_peddy(
    DATA_PATH,
    pseudobulk,
    dataset,
    ncells,
    rd,
    mod1="bulk_chunk_ribo",
    mod2="bulk_diss_polyA",
):
    """
    Read Peddy output and create a matrix for heatmap visualization.

    input:
        - DATA_PATH: base path to the data directory
        - pseudobulk: boolean indicating if the data is pseudobulk
        - dataset: the dataset name
        - ncells: the number of cells used to create the pseudobulk (or "null" for real data)
        - rd: the read depth filter cut-off used for the analysis
        - mod1: the first modality
        - mod2: the second modality

    output:
        - df: long format dataframe with columns [sample_a, sample_b, rel_difference]
        - matrix: square dataframe with samples as rows and columns, values are the relative differences between samples
    """
    if pseudobulk:
        rd_path = f"{DATA_PATH}/2a_peddy/{dataset}/pseudobulk/ncells_{ncells}/read_depth_{rd}"
    else:
        rd_path = _resolve_read_depth_path(
            DATA_PATH, "2a_peddy", dataset, pseudobulk, ncells,
            f"{mod1}_0_{mod2}_0", mod1, mod2,
        )
    file_path = os.path.join(rd_path, "output.ped_check.csv")
    df = pd.read_csv(
        file_path,
        header=0,
    )  
    df["sample_a"] = df["sample_a"].str.replace(".bam", "")
    df["sample_b"] = df["sample_b"].str.replace(".bam", "")

    df = df[["sample_a", "sample_b", "rel_difference"]]

    # peddy leaves out the diagonal values (self-comparisons) in the output, so we need to fill them in with 0s.
    samples = set(df["sample_a"].unique()) | set(df["sample_b"].unique())
    for sample in samples:
        if not df[(df["sample_a"] == sample) & (df["sample_b"] == sample)].empty:
            continue
        df = pd.concat([df, pd.DataFrame([{"sample_a": sample, "sample_b": sample, "rel_difference": 0}])], ignore_index=True)

    matrix = long_df_to_matrix_peddy(df)
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

    # Order the matrix rows and columns
    matrix = matrix.loc[sorted(matrix.index), sorted(matrix.columns)]

    return df, matrix


def parse_heatmap_matrix_somalier(
    DATA_PATH,
    pseudobulk,
    dataset,
    ncells,
    mod1="bulk_chunk_ribo",
    mod2="bulk_diss_polyA",
):
    """
    Read Somalier output and create a matrix for heatmap visualization.

    input:
        - DATA_PATH: base path to the data directory
        - pseudobulk: boolean indicating if the data is pseudobulk
        - dataset: the dataset name
        - ncells: the number of cells used to create the pseudobulk (or "null" for real data)
        - rd: the read depth filter cut-off used for the analysis
        - mod1: the first modality
        - mod2: the second modality

    output:
        - df: long format dataframe with columns [sample_a, sample_b, concordance]
        - matrix: square dataframe with samples as rows and columns, values are the concordance scores between samples
    """

    if pseudobulk:
        file_path = os.path.join(
            DATA_PATH,
            f"2b_somalier/{dataset}/pseudobulk/ncells_{ncells}/somalier/somalier.pairs.tsv",
        )
    else:
        base_path = _resolve_analysis_base_path(
            DATA_PATH,
            "2b_somalier",
            dataset,
            pseudobulk,
            ncells,
            mod1,
            mod2,
        )
        nested_file_path = os.path.join(base_path, "somalier", "somalier.pairs.tsv")
        flat_file_path = os.path.join(base_path, "somalier.pairs.tsv")
        file_path = (
            nested_file_path if os.path.exists(nested_file_path) else flat_file_path
        )
    df = pd.read_csv(
        file_path,
        sep="\t",
        header=0,
    )

    df["sample_a"] = df["#sample_a"].str.replace("_filtered", "")
    df["sample_b"] = df["sample_b"].str.replace("_filtered", "")
    df = df[["sample_a", "sample_b", "concordance"]]

    # We need to symmetrize the matrix by adding the reverse pairs, since Somalier only outputs one direction of the pairwise comparisons.
    df = pd.concat([df, df.rename(columns={"sample_a": "sample_b", "sample_b": "sample_a"})], ignore_index=True)
    # Somalier leaves out the diagonal values (self-comparisons) in the output, so we need to fill them in with 1s.
    samples = set(df["sample_a"].unique()) | set(df["sample_b"].unique())
    for sample in samples:
        if not df[(df["sample_a"] == sample) & (df["sample_b"] == sample)].empty:
            continue
        df = pd.concat([df, pd.DataFrame([{"sample_a": sample, "sample_b": sample, "concordance": 1}])], ignore_index=True)

    matrix = long_df_to_matrix_somalier(df, metric="concordance")

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

    # Order the matrix rows and columns
    matrix = matrix.loc[sorted(matrix.index), sorted(matrix.columns)]

    return df, matrix


def parse_heatmap_matrix_timeattackgencomp(
    DATA_PATH,
    pseudobulk,
    dataset,
    ncells,
    rd,
    mod1="bulk_chunk_ribo",
    mod2="bulk_diss_polyA",
):
    """
    Read TimeAttackGenComp output and create a matrix for heatmap visualization.

    input:
        - DATA_PATH: base path to the data directory
        - pseudobulk: boolean indicating if the data is pseudobulk
        - dataset: the dataset name
        - ncells: the number of cells used to create the pseudobulk (or "null" for real data)
        - rd: the read depth filter cut-off used for the analysis
        - mod1: the first modality
        - mod2: the second modality

    output:
        - df: long format dataframe with columns [sample1, sample2, value]
        - matrix: square dataframe with samples as rows and columns, values are the similarity scores used by TimeAttackGenComp
    """

    if pseudobulk:
        rd_path = f"{DATA_PATH}/2a_timeattackgencomp/{dataset}/pseudobulk/ncells_{ncells}/read_depth_{rd}"
    else:
        rd_path = _resolve_read_depth_path(
            DATA_PATH, "2a_timeattackgencomp", dataset, pseudobulk, ncells,
            f"{mod1}_0_{mod2}_0", mod1, mod2,
        )
    file_path = os.path.join(rd_path, "timeattackgencomp.snv.out.txt")

    matrix = pd.read_csv(file_path, sep="\t", header=0, index_col=0)[:-1]
    matrix.index = matrix.index.str.replace("_individual_variants", "")
    matrix.columns = matrix.columns.str.replace("_individual_variants", "")

    # Symmetrize matrix by copying the upper triangle to the lower triangle
    matrix = matrix.where(np.triu(np.ones(matrix.shape), k=0).astype(bool))
    matrix = matrix.combine_first(matrix.T)

    df = matrix.stack().reset_index()
    df.columns = ["sample1", "sample2", "value"]

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
    
    return df, matrix


def parse_heatmap_matrix_vireo(
    DATA_PATH,
    pseudobulk,
    dataset,
    ncells,
    rd,
    mod1="bulk_chunk_ribo",
    mod2="bulk_diss_polyA",
):
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
        rd_path = f"{DATA_PATH}/2a_vireo/{dataset}/pseudobulk/ncells_{ncells}/read_depth_{rd}/"
    else:
        rd_path = _resolve_read_depth_path(
            DATA_PATH, "2a_vireo", dataset, pseudobulk, ncells,
            f"{mod1}_0_{mod2}_0", mod1, mod2,
        )
    file_path = os.path.join(rd_path, "similarity_matrix.csv")
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

    if dataset == "hgsoc":
        matrix.columns = matrix.columns.str.replace(f"_{mod1}", "", 1).str.replace(f"_{mod2}", "", 1)
        matrix.index = matrix.index.str.replace(f"_{mod1}", "", 1).str.replace(f"_{mod2}", "", 1)

    # Vireo does not always return rows/columns for all samples, so we may have missing entries in the matrix
    # We add missing rows and columns with NaN values to ensure all samples are represented
    # Missing samples are identified using CrosscheckFingerprints
    crosscheck_df, _ = parse_heatmap_matrix_crosscheckfingerprints(
        DATA_PATH,
        pseudobulk,
        dataset,
        ncells,
        rd=0,
        mod1=mod1,
        mod2=mod2,
    )

    crosscheck_pairs = set()
    for _, row in crosscheck_df.iterrows():
        left, right = row["LEFT_SAMPLE"], row["RIGHT_SAMPLE"]
        crosscheck_pairs.add(tuple(sorted([left, right])))
    existing_pairs = set(
        tuple(sorted([i, j])) for i in matrix.index for j in matrix.columns
    )
    missing_pairs = crosscheck_pairs - existing_pairs
    for pair in missing_pairs:
        matrix.loc[pair[0], pair[1]] = np.nan
        matrix.loc[pair[1], pair[0]] = np.nan

    # Order rows and columns
    ordered_columns = sorted(matrix.columns)
    ordered_index = sorted(matrix.index)
    matrix = matrix.loc[ordered_index, ordered_columns]

    return matrix


#######################################################
### Functions to parse tool outputs and create sample match dataframes
#######################################################
def _load_standardized_heatmap_matrix(
    tool,
    dataset,
    ncells,
    read_depth,
    mod1,
    mod2,
    pseudobulk,
    standardized_heatmaps_path=None,
):
    """Load a numeric heatmap matrix previously exported by standardize_heatmap_data.py."""
    if standardized_heatmaps_path is None:
        standardized_heatmaps_path = os.path.abspath(
            os.path.join(
                os.path.dirname(__file__),
                "..",
                "benchmarking",
                "data",
                "standardized_heatmaps",
            )
        )

    if pseudobulk:
        if ncells is None:
            raise ValueError("ncells is required for pseudobulk heatmap loading")
        relative_path = os.path.join(dataset, "pseudobulk", f"ncells_{ncells}")
    else:
        relative_path = os.path.join(
            dataset, f"{mod1}_vs_{mod2}", "ncells_null"
        )

    if tool not in {"BAMixChecker", "Conpair", "Somalier", "ntsm"}:
        read_depth_component = str(read_depth).removeprefix("read_depth_")
        relative_path = os.path.join(
            relative_path, f"read_depth_{read_depth_component}"
        )

    file_path = os.path.join(standardized_heatmaps_path, relative_path, f"{tool}.csv")
    if not os.path.isfile(file_path):
        raise FileNotFoundError(
            f"Standardized heatmap not found for {tool}: {file_path}. "
            "Run standardize_heatmap_data.py first."
        )
    return pd.read_csv(file_path, index_col=0)


def parse_sample_matching_results_bamixchecker(
    DATA_PATH,
    pseudobulk,
    dataset,
    ncells,
    mod1="bulk_chunk_ribo",
    mod2="bulk_diss_polyA",
):
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
    df, _ = parse_heatmap_matrix_bamixchecker(
        DATA_PATH, pseudobulk, dataset, ncells, mod1, mod2
    )
    sample_matches = long_df_to_matrix_bamixchecker(
        df, metric="Conclusion"
    ).apply(lambda column: column.map(result_mapping))
    return sample_matches


def parse_sample_matching_results_conpair(
    DATA_PATH,
    pseudobulk,
    dataset,
    ncells,
    mod1="bulk_chunk_ribo",
    mod2="bulk_diss_polyA",
    threshold=0.8,
    standardized_heatmaps_path=None,
):
    """
    Parse Conpair output to determine sample match predictions.

    input: 
        - DATA_PATH: base path to the data directory
        - pseudobulk: boolean indicating if the data is pseudobulk
        - dataset: the dataset name
        - ncells: the number of cells used to create the pseudobulk (or "null" for real data)
        - mod1: the first modality
        - mod2: the second modality
        - threshold: concordance threshold for determining matches (default 0.8)

    output:
        - df: square dataframe with binary values indicating sample matches (1 for match, 0 for no match)
    """

    matrix = _load_standardized_heatmap_matrix(
        "Conpair",
        dataset,
        ncells,
        "0",
        mod1,
        mod2,
        pseudobulk,
        standardized_heatmaps_path,
    )
    sample_matches = (
        matrix.astype(float)
        .ge(threshold * 100)
        .where(matrix.notna())
        .astype("Int64")
    )

    return sample_matches


def parse_sample_matching_results_crosscheckfingerprints(
    DATA_PATH,
    pseudobulk,
    dataset,
    ncells,
    rd,
    mod1="bulk_chunk_ribo",
    mod2="bulk_diss_polyA",
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
        DATA_PATH, pseudobulk, dataset, ncells, rd, mod1, mod2
    )
    sample_matches = long_df_to_matrix_crosscheckfingerprints(df, metric="RESULT")
    sample_matches = sample_matches.apply(
        lambda column: column.map(result_mapping)
    )

    return sample_matches


def parse_sample_matching_results_hysys(
    DATA_PATH,
    pseudobulk,
    dataset,
    ncells,
    rd,
    mod1="bulk_chunk_ribo",
    mod2="bulk_diss_polyA",
):
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
    rd_path = _resolve_read_depth_path(
        DATA_PATH,
        "2a_hysys",
        dataset,
        pseudobulk,
        ncells,
        rd,
        mod1,
        mod2,
    )
    file_path = os.path.join(rd_path, "model_results.txt")

    with open(file_path, "r") as file:
        content = file.read()

    # Initialize collections
    related_pairs = set()
    inconclusive_pairs = set()
    all_samples = set()

    # Clean sample names function
    def clean_sample_name(sample_path):
        # Extract filename from path
        sample_name = sample_path.split("/")[-1]
        # Remove file extension
        sample_name = sample_name.replace(".snps", "")
        # Match parse_heatmap_matrix_* behavior: only strip dataset prefixes
        # for pseudobulk and hgsoc-new naming conventions.
        if pseudobulk or dataset == "hgsoc-new":
            regex_exp = dataset_regex_dict.get(dataset, "")
        else:
            regex_exp = ""
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
                related_pairs.add(pair)

    # Create all possible pairs from all samples
    all_pairs = [(s1, s2) for s1 in all_samples for s2 in all_samples]

    # Create DataFrame with all pairs
    pair_data = []
    for pair in all_pairs:
        if pair in inconclusive_pairs:
            match_status = float("nan")  # inconclusive
        elif pair in related_pairs or pair[::-1] in related_pairs:
            match_status = 1  # related pairs are predicted matches
        else:
            match_status = 0  # omitted pairs are predicted non-matches

        pair_data.append(
            {"sample1": pair[0], "sample2": pair[1], "match": match_status}
        )

    # Create matrix with sample1 as index and sample2 as columns
    df = pd.DataFrame(pair_data)
    matrix = df.pivot(index="sample1", columns="sample2", values="match")

    return matrix


def parse_sample_matching_results_ngscheckmate(
    DATA_PATH,
    pseudobulk,
    dataset,
    ncells,
    rd,
    mod1="bulk_chunk_ribo",
    mod2="bulk_diss_polyA",
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
        DATA_PATH, pseudobulk, dataset, ncells, rd, mod1, mod2
    )
    sample_matches = long_df_to_matrix_ngscheckmate(
        df, "Matched"
    ).apply(lambda column: column.map(result_mapping))

    return sample_matches


def parse_sample_matching_results_ntsm(
    DATA_PATH,
    pseudobulk,
    dataset,
    ncells,
    mod1="bulk_chunk_ribo",
    mod2="bulk_diss_polyA",
    standardized_heatmaps_path=None,
):
    """
    Parse NTSM results to categorize sample relationships.

    input:
        - DATA_PATH: base path to the data directory
        - pseudobulk: boolean indicating if the data is pseudobulk
        - dataset: the dataset name
        - ncells: the number of cells used to create the pseudobulk (or "null" for real data)

    output:
        - df: square dataframe with binary values indicating sample matches (1 for match, 0 for no match)
    """
    matrix = _load_standardized_heatmap_matrix(
        "ntsm",
        dataset,
        ncells,
        "0",
        mod1,
        mod2,
        pseudobulk,
        standardized_heatmaps_path,
    )
    # If the ntsm relatedness is greater than zero, then samples are considered a match (1), otherwise not a match (0)
    sample_matches = (matrix > 0).astype(int)

    return sample_matches


def parse_sample_matching_results_omicsprint(
    DATA_PATH,
    pseudobulk,
    dataset,
    ncells,
    rd,
    mod1="bulk_chunk_ribo",
    mod2="bulk_diss_polyA",
    threshold=1.7,
    standardized_heatmaps_path=None,
):
    """
    Parse OmicsPrint output to determine sample match predictions.

    input: 
        - DATA_PATH: base path to the data directory
        - pseudobulk: boolean indicating if the data is pseudobulk
        - dataset: the dataset name
        - ncells: the number of cells used to create the pseudobulk (or "null" for real data)
        - rd: the read depth filter cut-off used for the analysis
        - mod1: the first modality
        - mod2: the second modality
        - threshold: concordance threshold for determining matches (default 1.7)

    output:
        - df: square dataframe with binary values indicating sample matches (1 for match, 0 for no match)
    """
    matrix = _load_standardized_heatmap_matrix(
        "OmicsPrint",
        dataset,
        ncells,
        rd,
        mod1,
        mod2,
        pseudobulk,
        standardized_heatmaps_path,
    )
    # If the omicsprint value is >= threshold, then samples are considered a match (1), otherwise not a match (0).
    # NaNs are preserved
    sample_matches = (
        matrix.ge(threshold)
        .where(matrix.notna())
        .astype("Int64")
    )

    return sample_matches


def parse_sample_matching_results_peddy(
    DATA_PATH,
    pseudobulk,
    dataset,
    ncells,
    rd,
    mod1="bulk_chunk_ribo",
    mod2="bulk_diss_polyA",
    threshold=0,
    standardized_heatmaps_path=None,
):
    """
    Parse Peddy output to categorize sample relationships.
    
    input:
        - DATA_PATH: base path to the data directory
        - pseudobulk: boolean indicating if the data is pseudobulk
        - dataset: the dataset name
        - ncells: the number of cells used to create the pseudobulk (or "null" for real data)
        - rd: the read depth filter cut-off used for the analysis
        - mod1: the first modality
        - mod2: the second modality
        - threshold: the relative difference threshold for determining matches (default 0)

    output:
        - df: square dataframe with binary values indicating sample matches (1 for match, 0 for no match)
    """

    matrix = _load_standardized_heatmap_matrix(
        "Peddy",
        dataset,
        ncells,
        rd,
        mod1,
        mod2,
        pseudobulk,
        standardized_heatmaps_path,
    )
    # If peddy values are <= threshold, then samples are considered a match (1), otherwise not a match (0)
    # NaNs are preserved
    sample_matches = (
        matrix.le(threshold)
        .where(matrix.notna())
        .astype("Int64")
    )
    
    return sample_matches


def parse_sample_matching_results_somalier(
    DATA_PATH,
    pseudobulk,
    dataset,
    ncells,
    mod1="bulk_chunk_ribo",
    mod2="bulk_diss_polyA",
    threshold=1,
    standardized_heatmaps_path=None,
):
    """
    Parse Somalier output to categorize sample relationships.

    input:
        - DATA_PATH: base path to the data directory
        - pseudobulk: boolean indicating if the data is pseudobulk
        - dataset: the dataset name
        - ncells: the number of cells used to create the pseudobulk (or "null" for real data)
        - mod1: the first modality
        - mod2: the second modality
        - threshold: the concordance threshold for determining matches (default 0)

    output:
        - df: square dataframe with binary values indicating sample matches (1 for match, 0 for no match)
    """
    matrix = _load_standardized_heatmap_matrix(
        "Somalier",
        dataset,
        ncells,
        "0",
        mod1,
        mod2,
        pseudobulk,
        standardized_heatmaps_path,
    )
    # If somalier value is >= threshold, then samples are considered a match (1), otherwise not a match (0)
    # NaNs are preserved
    sample_matches = (
        matrix.ge(threshold)
        .where(matrix.notna())
        .astype("Int64")
    )

    return sample_matches


def parse_sample_matching_results_timeattackgencomp(
    DATA_PATH,
    pseudobulk,
    dataset,
    ncells,
    read_depth,
    mod1="bulk_chunk_ribo",
    mod2="bulk_diss_polyA",
    threshold=0.05,
    standardized_heatmaps_path=None,
):  
    """
    Parse TimeAttackGenComp output to categorize sample relationships.

    input:
        - DATA_PATH: base path to the data directory
        - pseudobulk: boolean indicating if the data is pseudobulk
        - dataset: the dataset name
        - ncells: the number of cells used to create the pseudobulk (or "null" for real data)
        - read_depth: the read depth filter cut-off used for the analysis
        - mod1: the first modality
        - mod2: the second modality

    output:
        - df: square dataframe with binary values indicating sample matches (1 for match, 0 for no match)
    """

    matrix = _load_standardized_heatmap_matrix(
        "TimeAttackGenComp",
        dataset,
        ncells,
        read_depth,
        mod1,
        mod2,
        pseudobulk,
        standardized_heatmaps_path,
    )
    # If timeattackgencomp value is 0, then samples are considered a match (1), otherwise not a match (0)
    # NaNs are preserved
    sample_matches = (
        matrix.le(threshold)
        .where(matrix.notna())
        .astype("Int64")
    )

    return sample_matches


def parse_sample_matching_results_vireo(
    DATA_PATH,
    pseudobulk,
    dataset,
    ncells,
    rd,
    mod1="bulk_chunk_ribo",
    mod2="bulk_diss_polyA",
    standardized_heatmaps_path=None,
):
    """
    Parse Vireo predictions using the standardized heatmap sample universe.

    Real-data matches come from matched_samples.csv. Pseudobulk matches are
    inferred by minimizing the standardized Vireo similarity matrix with the
    linear sum assignment algorithm.

    output:
        - rectangular dataframe with 1 for matches, 0 for non-matches, and NaN
          wherever the standardized heatmap has an unknown score
    """
    heatmap = _load_standardized_heatmap_matrix(
        "Vireo",
        dataset,
        ncells,
        rd,
        mod1,
        mod2,
        pseudobulk,
        standardized_heatmaps_path,
    ).astype(float)

    predictions = pd.DataFrame(
        0.0, index=heatmap.index.copy(), columns=heatmap.columns.copy()
    )

    if pseudobulk:
        predictions = predictions.where(heatmap.notna())
        scores = heatmap.to_numpy(dtype=float)
        finite = np.isfinite(scores)
        if finite.any():
            penalty = max(1.0, float(np.max(np.abs(scores[finite])))) * 1e9
            assignment_cost = np.where(finite, scores, penalty)
            row_indices, column_indices = linear_sum_assignment(assignment_cost)
            for row_index, column_index in zip(row_indices, column_indices):
                if finite[row_index, column_index]:
                    predictions.iat[row_index, column_index] = 1.0
        return predictions

    rd_path = _resolve_read_depth_path(
        DATA_PATH,
        "2a_vireo",
        dataset,
        pseudobulk,
        ncells,
        rd,
        mod1,
        mod2,
    )
    file_path = os.path.join(rd_path, "matched_samples.csv")
    matched_pairs = pd.read_csv(file_path)
    regex_exp = dataset_regex_dict.get(dataset, "") if dataset == "hgsoc-new" else ""

    def normalize_sample_label(sample):
        normalized = re.sub(
            regex_exp,
            "",
            str(sample).replace(".bam", ""),
        )
        normalized = normalized.replace("ds.", "")
        return normalized.replace("-", "_")

    matched_pairs = {
        (
            normalize_sample_label(sample0),
            normalize_sample_label(sample1),
        )
        for sample0, sample1 in matched_pairs[["sample_id_0", "sample_id_1"]].itertuples(
            index=False, name=None
        )
    }
    normalized_columns = {
        column: normalize_sample_label(column) for column in predictions.columns
    }
    for row_label in predictions.index:
        normalized_row = normalize_sample_label(row_label)
        for column_label, normalized_column in normalized_columns.items():
            pair = (normalized_row, normalized_column)
            reverse_pair = (normalized_column, normalized_row)
            if pair in matched_pairs or reverse_pair in matched_pairs:
                predictions.loc[row_label, column_label] = 1.0

    return predictions


def load_sample_matching_results(
    DATA_PATH,
    pseudobulk,
    tool,
    dataset,
    ncells,
    rd1,
    rd2=None,
    mod1="bulk_chunk_ribo",
    mod2="bulk_diss_polyA",
    threshold=None,
    standardized_heatmaps_path=None,
):

    if rd2 is not None:
        rd = f"{mod1}_{rd1}_{mod2}_{rd2}"
        read_depth_label = f"{mod1}={rd1}, {mod2}={rd2}"
    else:
        rd = rd1
        read_depth_label = str(rd1)

    print(
        f"Processing tool {tool}, dataset {dataset}, read depth {read_depth_label}..."
    )
    try:
        if tool == "BAMixChecker":
            sample_matches = parse_sample_matching_results_bamixchecker(
                DATA_PATH, pseudobulk, dataset, ncells, mod1, mod2
            )
        elif tool == "Conpair":
            parser_args = (
                DATA_PATH,
                pseudobulk,
                dataset,
                ncells,
                mod1,
                mod2,
            )
            sample_matches = (
                parse_sample_matching_results_conpair(
                    *parser_args,
                    standardized_heatmaps_path=standardized_heatmaps_path,
                )
                if threshold is None
                else parse_sample_matching_results_conpair(
                    *parser_args,
                    threshold=threshold,
                    standardized_heatmaps_path=standardized_heatmaps_path,
                )
            )
        elif tool == "CrosscheckFingerprints":
            sample_matches = parse_sample_matching_results_crosscheckfingerprints(
                DATA_PATH, pseudobulk, dataset, ncells, rd, mod1, mod2
            )
        elif tool == "HYSYS":
            sample_matches = parse_sample_matching_results_hysys(
                DATA_PATH, pseudobulk, dataset, ncells, rd, mod1, mod2
            )
        elif tool == "NGSCheckmate":
            sample_matches = parse_sample_matching_results_ngscheckmate(
                DATA_PATH, pseudobulk, dataset, ncells, rd, mod1, mod2
            )
        elif tool == "ntsm":
            sample_matches = parse_sample_matching_results_ntsm(
                DATA_PATH,
                pseudobulk,
                dataset,
                ncells,
                mod1,
                mod2,
                standardized_heatmaps_path=standardized_heatmaps_path,
            )
        elif tool == "OmicsPrint":
            parser_args = (DATA_PATH, pseudobulk, dataset, ncells, rd, mod1, mod2)
            sample_matches = (
                parse_sample_matching_results_omicsprint(
                    *parser_args,
                    standardized_heatmaps_path=standardized_heatmaps_path,
                )
                if threshold is None
                else parse_sample_matching_results_omicsprint(
                    *parser_args,
                    threshold=threshold,
                    standardized_heatmaps_path=standardized_heatmaps_path,
                )
            )
        elif tool == "Peddy":
            parser_args = (DATA_PATH, pseudobulk, dataset, ncells, rd, mod1, mod2)
            sample_matches = (
                parse_sample_matching_results_peddy(
                    *parser_args,
                    standardized_heatmaps_path=standardized_heatmaps_path,
                )
                if threshold is None
                else parse_sample_matching_results_peddy(
                    *parser_args,
                    threshold=threshold,
                    standardized_heatmaps_path=standardized_heatmaps_path,
                )
            )
        elif tool == "Somalier":
            parser_args = (DATA_PATH, pseudobulk, dataset, ncells, mod1, mod2)
            sample_matches = (
                parse_sample_matching_results_somalier(
                    *parser_args,
                    standardized_heatmaps_path=standardized_heatmaps_path,
                )
                if threshold is None
                else parse_sample_matching_results_somalier(
                    *parser_args,
                    threshold=threshold,
                    standardized_heatmaps_path=standardized_heatmaps_path,
                )
            )
        elif tool == "TimeAttackGenComp":
            parser_args = (DATA_PATH, pseudobulk, dataset, ncells, rd, mod1, mod2)
            sample_matches = (
                parse_sample_matching_results_timeattackgencomp(
                    *parser_args,
                    standardized_heatmaps_path=standardized_heatmaps_path,
                )
                if threshold is None
                else parse_sample_matching_results_timeattackgencomp(
                    *parser_args,
                    threshold=threshold,
                    standardized_heatmaps_path=standardized_heatmaps_path,
                )
            )
        elif tool == "Vireo":
            sample_matches = parse_sample_matching_results_vireo(
                DATA_PATH,
                pseudobulk,
                dataset,
                ncells,
                rd,
                mod1,
                mod2,
                standardized_heatmaps_path=standardized_heatmaps_path,
            )
        else:
            print(f"Error! Do not recognize tool {tool}")
            print(
                "Choice of tool must be one of BAMixChecker, Conpair, CrosscheckFingerprints, HYSYS, NGSCheckmate, ntsm,OmicsPrint, Peddy, Somalier, TimeAttackGenComp, or Vireo."
            )
            sample_matches = None
    except FileNotFoundError:
        print(
            f"Could not find data for {tool}, real dataset {dataset}, read depth {read_depth_label}. "
        )
        sample_matches = None

    if tool != "Vireo" and sample_matches is not None:
        # Filter sample_matches to only include mod1 samples in the first column and mod2 samples in the second column
        sample_matches = sample_matches.loc[
            [idx for idx in sample_matches.index if _modality_in_label(mod1, idx)],
            [col for col in sample_matches.columns if _modality_in_label(mod2, col)],
        ]

    return sample_matches


def load_heatmap_data(
    DATA_PATH,
    pseudobulk,
    tool,
    dataset,
    ncells,
    rd,
    mod1="bulk_chunk_ribo",
    mod2="bulk_diss_polyA",
):

    print(f"Processing tool {tool}, dataset {dataset}, read depth {rd}...")
    try:
        # Create heatmap matrix
        if tool == "BAMixChecker":
            _, matrix = parse_heatmap_matrix_bamixchecker(
                DATA_PATH, pseudobulk, dataset, ncells, mod1, mod2
            )
        elif tool == "Conpair":
            _, matrix = parse_heatmap_matrix_conpair(
                DATA_PATH, pseudobulk, dataset, ncells, mod1, mod2
            )
        elif tool == "CrosscheckFingerprints":
            _, matrix = parse_heatmap_matrix_crosscheckfingerprints(
                DATA_PATH, pseudobulk, dataset, ncells, rd, mod1, mod2
            )
        elif tool == "HYSYS":
            _, matrix = parse_heatmap_matrix_hysys(
                DATA_PATH, pseudobulk, dataset, ncells, rd, mod1, mod2
            )
        elif tool == "NGSCheckmate":
            _, matrix = parse_heatmap_matrix_ngscheckmate(
                DATA_PATH, pseudobulk, dataset, ncells, rd, mod1, mod2
            )
        elif tool == "ntsm":
            _, matrix = parse_heatmap_matrix_ntsm(
                DATA_PATH, pseudobulk, dataset, ncells, mod1, mod2
            )
        elif tool == "OmicsPrint":
            _, matrix = parse_heatmap_matrix_omicsprint(
                DATA_PATH, pseudobulk, dataset, ncells, rd, mod1, mod2
            )
        elif tool == "Peddy":
            _, matrix = parse_heatmap_matrix_peddy(
                DATA_PATH, pseudobulk, dataset, ncells, rd, mod1, mod2
            )
        elif tool == "Somalier":
            _, matrix = parse_heatmap_matrix_somalier(
                DATA_PATH, pseudobulk, dataset, ncells, mod1, mod2
            )
        elif tool == "TimeAttackGenComp":
            _, matrix = parse_heatmap_matrix_timeattackgencomp(
                DATA_PATH, pseudobulk, dataset, ncells, rd, mod1, mod2
            )
        elif tool == "Vireo":
            matrix = parse_heatmap_matrix_vireo(
                DATA_PATH, pseudobulk, dataset, ncells, rd, mod1, mod2
            )
        else:
            print(f"Error! Do not recognize tool {tool}")
            print(
                "Choice of tool must be one of BAMixChecker, Conpair, CrosscheckFingerprints, HYSYS, NGSCheckmate, ntsm, OmicsPrint, Peddy, Somalier, TimeAttackGenComp, or Vireo."
            )
            matrix = None
    except FileNotFoundError:
        print(
            f"Could not find data for {tool}, real dataset {dataset}, read depth {rd}. "
        )
        matrix = None

    return matrix


def order_matrix_by_expected_matches(
    matrix, expected_matches, mod1="bulk", mod2="single-cell"
):
    """
    Add any missing samples from mod1 and mod2 to the matrix as NaN rows/columns, then
    order the rows and columns of the matrix according to the expected matches between mod1 and mod2 samples,
    so that samples that are expected to match show up on the diagonal of the heatmap.
    Samples with no expected match will be shown after the expected matches

    input:
        - matrix: dataframe containing similarity measures of mod1 samples (rows) against mod2 samples (columns)
        - expected_matches: dataframe containing expected matches between mod1 and mod2 samples, with one column for mod1 sample names and one column for mod2 sample names
        - mod1: name of first modality (e.g. "bulk")
        - mod2: name of second modality (e.g. "single-cell")

    output:
        - matrix ordered according to expected matches
    """
    def _norm(text):
        return str(text).lower().replace("-", " ").replace("_", " ").strip()

    def _resolve_expected_modality_columns_local(matches_df, modality1, modality2):
        def _has_single_family(text):
            key = _norm(text)
            return any(token in key for token in ["single", "nucleus", "pooled"])

        columns = list(matches_df.columns)
        norm_cols = {_norm(col): col for col in columns}

        mod1_key = _norm(modality1)
        mod2_key = _norm(modality2)

        mod1_candidates = [orig for norm, orig in norm_cols.items() if mod1_key in norm]
        mod2_candidates = [orig for norm, orig in norm_cols.items() if mod2_key in norm]

        mod1_col = mod1_candidates[0] if len(mod1_candidates) >= 1 else None
        mod2_col = mod2_candidates[0] if len(mod2_candidates) >= 1 else None

        # Conservative fallback: only infer the counterpart column when the unresolved
        # requested modality is in the single/pooled family and the remaining column is too.
        if mod1_col is not None and mod2_col is None and len(columns) == 2:
            other_col = [col for col in columns if col != mod1_col][0]
            if _has_single_family(modality2) and _has_single_family(other_col):
                mod2_col = other_col
        if mod2_col is not None and mod1_col is None and len(columns) == 2:
            other_col = [col for col in columns if col != mod2_col][0]
            if _has_single_family(modality1) and _has_single_family(other_col):
                mod1_col = other_col

        if mod1_col is None or mod2_col is None:
            raise ValueError(
                "Could not resolve expected match columns for requested modalities "
                f"mod1={modality1}, mod2={modality2}. Available columns: {columns}"
            )

        return mod1_col, mod2_col

    mod1_col, mod2_col = _resolve_expected_modality_columns_local(
        expected_matches, mod1, mod2
    )

    mod1_series = expected_matches[mod1_col]
    mod2_series = expected_matches[mod2_col]

    # Some result files store matrix axes in the opposite orientation (mod2 rows, mod1 columns).
    # Detect and correct this before attempting expected-match alignment.
    row_has_mod1 = sum(_modality_in_label(mod1, label) for label in matrix.index)
    row_has_mod2 = sum(_modality_in_label(mod2, label) for label in matrix.index)
    col_has_mod1 = sum(_modality_in_label(mod1, label) for label in matrix.columns)
    col_has_mod2 = sum(_modality_in_label(mod2, label) for label in matrix.columns)
    if (
        row_has_mod1 == 0
        and col_has_mod2 == 0
        and col_has_mod1 > 0
        and row_has_mod2 > 0
    ):
        matrix = matrix.T

    # Keep exact expected ordering; make NaN positions explicit and unique instead of collapsing to "X".
    mod1_samples = [
        str(value) if pd.notna(value) else f"__{mod1}_missing_expected_{i}"
        for i, value in enumerate(mod1_series)
    ]
    mod2_samples = [
        str(value) if pd.notna(value) else f"__{mod2}_missing_expected_{i}"
        for i, value in enumerate(mod2_series)
    ]

    def sample_matches_label(sample, label):
        sample_str = re.escape(str(sample))
        return bool(
            re.search(
                rf"(^|[_/\\-]){sample_str}($|[_/\\-])",
                str(label),
                flags=re.IGNORECASE,
            )
        )

    # Match expected rows/cols to source matrix labels one-by-one, consuming each source label at most once.
    # Restrict matches to the correct modality to avoid cross-modality mismatches.
    available_rows = [
        label for label in matrix.index if _modality_in_label(mod1, label)
    ]
    matched_rows = []
    for sample in mod1_samples:
        row_match = None
        for idx, label in enumerate(available_rows):
            if sample_matches_label(sample, label):
                row_match = label
                del available_rows[idx]
                break
        matched_rows.append(row_match)

    available_cols = [
        label for label in matrix.columns if _modality_in_label(mod2, label)
    ]
    matched_cols = []
    for sample in mod2_samples:
        col_match = None
        for idx, label in enumerate(available_cols):
            if sample_matches_label(sample, label):
                col_match = label
                del available_cols[idx]
                break
        matched_cols.append(col_match)

    # Build output axis labels in expected order.
    # - matched entries: keep their source matrix label
    # - missing entries: "{sample}_missing" for known samples, "missing_{i}" for NaN entries
    output_row_labels = []
    for i, (sample, source_row) in enumerate(zip(mod1_samples, matched_rows)):
        if source_row is not None:
            label = source_row
        elif sample.startswith("__"):  # was NaN in expected_matches
            label = f"missing_{i}"
        else:
            label = f"{sample}_missing"
        output_row_labels.append(label)

    output_col_labels = []
    for i, (sample, source_col) in enumerate(zip(mod2_samples, matched_cols)):
        if source_col is not None:
            label = source_col
        elif sample.startswith("__"):  # was NaN in expected_matches
            label = f"missing_{i}"
        else:
            label = f"{sample}_missing"
        output_col_labels.append(label)

    # Start from an empty matrix with the exact expected ordering.
    ordered_matrix = pd.DataFrame(
        np.nan, index=output_row_labels, columns=output_col_labels
    )

    # Fill values where both expected row and expected column were found in the source matrix.
    for row_pos, source_row in enumerate(matched_rows):
        if source_row is None:
            continue
        for col_pos, source_col in enumerate(matched_cols):
            if source_col is None:
                continue
            ordered_matrix.iat[row_pos, col_pos] = matrix.loc[source_row, source_col]

    return ordered_matrix
