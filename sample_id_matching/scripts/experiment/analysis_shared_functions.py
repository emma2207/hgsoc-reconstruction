import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


dataset_regex_dict = {
    "hgsoc": r"[a-zA-Z0-9]{32}_",
    "low_grade_glioma": r"^GSM[0-9]{7}_"
}


def long_df_to_matrix_bamixchecker(df, metric="Concordance Rate"):
    # Use pivot_table for much faster matrix creation
    matrix = df.pivot_table(
        index=0,
        columns='Sample',
        values=metric,
        aggfunc='first'  # Use first value if duplicates exist
    )
    
    # Get all unique samples to create a square matrix
    all_samples = sorted(set(df.index) | set(df['Sample']))
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
        index='LEFT_SAMPLE',
        columns='RIGHT_SAMPLE', 
        values=metric,
        aggfunc='first'  # Use first value if duplicates exist
    )
    
    # Get all unique samples to create a square matrix
    all_samples = sorted(set(df['LEFT_SAMPLE']) | set(df['RIGHT_SAMPLE']))
    matrix = matrix.reindex(index=all_samples, columns=all_samples)
    
    # Make matrix symmetric by filling NaN values
    matrix = matrix.combine_first(matrix.T)
    
    # Order rows and columns by bulk modality and sample ID
    matrix = matrix.sort_index(
        key=lambda x: x.str.split('_').str[0]).sort_index(key=lambda x: x.str.split('_', n=1).str[1]
    )
    matrix = matrix.reindex(
        sorted(
            sorted(
                matrix.columns,
                key=lambda x: x.split('_', 1)[0]
            ),
            key=lambda x: x.split('_', 1)[1]
        )
    )
    matrix = matrix[matrix.index]

    return matrix


def long_df_to_matrix_hysys(df):
    # Use pivot_table for much faster matrix creation
    matrix = df.pivot_table(
        index=0,
        columns='Sample',
        values='Concordance',
        aggfunc='first'  # Use first value if duplicates exist
    )
    
    # Get all unique samples to create a square matrix
    all_samples = sorted(set(df.index) | set(df['Sample']))
    matrix = matrix.reindex(index=all_samples, columns=all_samples)
    
    # Make matrix symmetric by filling NaN values
    matrix = matrix.combine_first(matrix.T)
    
    # Order rows and columns by bulk modality and sample ID
    matrix = matrix.sort_index(
        key=lambda x: x.str.split('_').str[0]).sort_index(key=lambda x: x.str.split('_', n=1).str[1]
    )
    matrix = matrix.reindex(
        sorted(
            sorted(
                matrix.columns,
                key=lambda x: x.split('_', 1)[0]
            ),
            key=lambda x: x.split('_', 1)[1]
        )
    )
    matrix = matrix[matrix.index]

    return matrix


def long_df_to_matrix_ngscheckmate(df, metric="Correlation"):
    # Use pivot_table for much faster matrix creation
    matrix = df.pivot_table(
        index=0,
        columns='Sample',
        values=metric,
        aggfunc='first'  # Use first value if duplicates exist
    )
    
    # Convert correlation values to percentage (divide by 100)
    if metric == "Correlation":
        matrix = matrix / 100
    
    # Get all unique samples to create a square matrix
    all_samples = sorted(set(df.index) | set(df['Sample']))
    matrix = matrix.reindex(index=all_samples, columns=all_samples)
    
    # Make matrix symmetric by filling NaN values
    matrix = matrix.combine_first(matrix.T)
    
    # Order rows and columns by bulk modality and sample ID
    matrix = matrix.sort_index(
        key=lambda x: x.str.split('_').str[0]).sort_index(key=lambda x: x.str.split('_', n=1).str[1]
    )
    matrix = matrix.reindex(
        sorted(
            sorted(
                matrix.columns,
                key=lambda x: x.split('_', 1)[0]
            ),
            key=lambda x: x.split('_', 1)[1]
        )
    )
    matrix = matrix[matrix.index]

    return matrix


def parse_heatmap_matrix_bamixchecker(DATA_PATH, pseudobulk, dataset, ncells):
    if pseudobulk:
        df = pd.read_csv(
            DATA_PATH + 
            f"2b_bamixchecker/{dataset}/pseudobulk/ncells_{ncells}/BAMixChecker/Total_result.txt", 
            sep="\t", index_col=0, header=None
        )
    else:
        df = pd.read_csv(
            DATA_PATH + 
            f"2b_bamixchecker/{dataset}/real_data/ncells_null/BAMixChecker/Total_result.txt", 
            sep="\t", index_col=0, header=None
        )
        
    df.columns = ["Sample", "Concordance Rate", "Conclusion"]

    regex_exp = dataset_regex_dict.get(dataset, "")
    rename_dict = {
        x: 
        re.sub(regex_exp, "", x.replace("_filtered_rg.bam", ""))
        for x in list(set(df.index) | set(df["Sample"]))
    }
    df.index = df.index.map(rename_dict)
    df["Sample"] = df["Sample"].map(rename_dict)
    
    matrix = long_df_to_matrix_bamixchecker(df)

    return df, matrix


def parse_heatmap_matrix_crosscheckfingerprints(DATA_PATH, pseudobulk, dataset, ncells, rd):
    # Read CrosscheckFingerprints output
    if pseudobulk:
        df = pd.read_csv(
            DATA_PATH + 
            f"2a_fingerprints/{dataset}/pseudobulk/ncells_{ncells}/read_depth_{rd}/crosscheck_metrics.txt", 
            sep="\t", skiprows=6, header=0
        )
    else:
        df = pd.read_csv(
            DATA_PATH + 
            f"2a_fingerprints/{dataset}/real_data/ncells_null/read_depth_{rd}/crosscheck_metrics.txt", 
            sep="\t", skiprows=6, header=0
        )
    df = df[["LEFT_SAMPLE", "RIGHT_SAMPLE", "LOD_SCORE", "RESULT"]]

    regex_exp = dataset_regex_dict.get(dataset, "")

    rename_dict = {
        x: 
        re.sub(regex_exp, "", x.replace(".bam", "").replace("ds.", ""))
        
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
            DATA_PATH + f"2a_hysys/{dataset}/pseudobulk/ncells_{ncells}/read_depth_{rd}/concordance_output.txt", 
            sep="\t", index_col=0, header=None
        )
    else:
        df = pd.read_csv(
            DATA_PATH + f"2a_hysys/{dataset}/real_data/ncells_null/read_depth_{rd}/concordance_output.txt", 
            sep="\t", index_col=0, header=None
        )

    regex_exp = dataset_regex_dict.get(dataset, "")

    rename_dict = {
        x: 
        re.sub(regex_exp, "", x.replace(f"{dataset}/", "").replace("_individual_variants.snps", "").replace("pseudobulk/", "").replace("real_data/", "").replace(f"ncells_{ncells}/", "").replace(f"read_depth_{rd}/", ""))
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
            DATA_PATH + 
            f"2a_ngscheckmate/{dataset}/pseudobulk/ncells_{ncells}/read_depth_{rd}/output_all.txt", 
            sep="\t", index_col=0, header=None
        )
    else:
        df = pd.read_csv(
            DATA_PATH + 
            f"2a_ngscheckmate/{dataset}/real_data/ncells_{ncells}/read_depth_{rd}/output_all.txt", 
            sep="\t", index_col=0, header=None
        )

    regex_exp = dataset_regex_dict.get(dataset, "")

    # Clean up sample names
    rename_dict = {
        x: 
        re.sub(regex_exp, "", x.replace(".vcf", "").replace("ds.", "").replace("_individual_variants", ""))
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
            DATA_PATH + 
            f"2a_vireo/{dataset}/pseudobulk/ncells_{ncells}/read_depth_{rd}/similarity_matrix.csv", 
            index_col=0
        )
    else:
        matrix = pd.read_csv(
            DATA_PATH + 
            f"2a_vireo/{dataset}/real_data/ncells_{ncells}/read_depth_{rd}/similarity_matrix.csv", 
            index_col=0
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


def parse_hysys_model_results(file_path, dataset):
    """
    Parse HYSYS model_results.txt file to categorize sample relationships.
    
    Args:
        file_path (str): Path to the model_results.txt file
        
    Returns:
        pd.DataFrame: DataFrame with MultiIndex (sample1, sample2) and 'match' column
                     Values: 1 (matching), 0 (not matching), NaN (inconclusive)
    """
    import re
    import pandas as pd
    import itertools
    
    with open(file_path, 'r') as file:
        content = file.read()
    
    # Initialize collections
    matching_pairs = set()
    inconclusive_pairs = set()
    all_samples = set()
    
    # Clean sample names function
    def clean_sample_name(sample_path):
        # Extract filename from path
        sample_name = sample_path.split('/')[-1]
        # Remove file extension
        sample_name = sample_name.replace('.snps', '')
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
        related_samples = [s.strip().strip("'\"") for s in related_samples_str.split(',')]
        
        for sample2_path in related_samples:
            sample2 = clean_sample_name(sample2_path.strip())
            all_samples.add(sample2)
            if sample1 != sample2:  # Don't include self-relationships
                pair = tuple(sorted([sample1, sample2]))  # Sort to avoid duplicates
                matching_pairs.add(pair)
    
    # Create all possible pairs from all samples
    all_samples_list = sorted(list(all_samples))
    all_pairs = list(itertools.combinations(all_samples_list, 2))
    
    # Create DataFrame with all pairs
    pair_data = []
    for pair in all_pairs:
        if pair in matching_pairs:
            match_status = 1  # matching
        elif pair in inconclusive_pairs:
            match_status = float('nan')  # inconclusive
        else:
            match_status = 0  # not matching
        
        pair_data.append({
            'sample1': pair[0],
            'sample2': pair[1], 
            'match': match_status
        })
    
    # Create DataFrame with MultiIndex
    df = pd.DataFrame(pair_data)
    df = df.set_index(['sample1', 'sample2'])
    
    return df

def print_hysys_results_summary(results_df):
    """Print a summary of HYSYS results."""
    matching_count = (results_df['match'] == 1).sum()
    not_matching_count = (results_df['match'] == 0).sum()
    inconclusive_count = results_df['match'].isna().sum()
    
    print("HYSYS Model Results Summary:")
    print("=" * 40)
    print(f"Matching pairs: {matching_count}")
    print(f"Not matching pairs: {not_matching_count}")
    print(f"Inconclusive pairs: {inconclusive_count}")
    print(f"Total pairs: {len(results_df)}")
    print()
    
    if matching_count > 0:
        matching_pairs = results_df[results_df['match'] == 1]
        print("Sample pairs that match:")
        for i, (idx, row) in enumerate(matching_pairs.head(10).iterrows()):
            print(f"  {idx[0]} ↔ {idx[1]}")
        if len(matching_pairs) > 10:
            print(f"  ... and {len(matching_pairs) - 10} more")
        print()
    
    if inconclusive_count > 0:
        inconclusive_pairs = results_df[results_df['match'].isna()]
        print("Sample pairs that are inconclusive:")
        for i, (idx, row) in enumerate(inconclusive_pairs.head(10).iterrows()):
            print(f"  {idx[0]} ↔ {idx[1]}")
        if len(inconclusive_pairs) > 10:
            print(f"  ... and {len(inconclusive_pairs) - 10} more")


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
    ax.set_title(f"{tool} Similarity Matrix", pad=20)

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
    sub_matrix = matrix_df[
        [
            col
            for col in matrix_df.columns
            if "single-cell" in col or "single-nucleus" in col
        ]
    ]
    sub_matrix = sub_matrix.loc[[idx for idx in sub_matrix.index if "bulk" in idx]]

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
    ax.set_ylabel("Bulk samples")
    ax.set_title(f"{tool} Bulk vs {sc} Similarity Matrix", pad=20)

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
