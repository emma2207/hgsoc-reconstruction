import os
import pandas as pd
from typing import Literal


def prepare_count_matrix(INPUT_PATH: str) -> pd.DataFrame:
    """From the folder with alignment results, this function traverses
    the folder hierarchy to find all the readspergene.out.tab files and
    combines them into a count matrix.
    """

    df = pd.DataFrame()

    for subdir, _, files in os.walk(INPUT_PATH):
        for file in files:
            if file.endswith("out.tab"):
                single_sample_counts = pd.read_csv(
                    subdir + "/" + file, delimiter="\t", header=[0, 1, 2, 3]
                )
                # Rename columns
                single_sample_counts.columns = [
                    "Gene ID",
                    "Counts Unstranded RNA-seq",
                    "Counts 1st Read Strand Aligned with RNA",
                    "Counts 2nd Read Strand Aligned with RNA",
                ]
                # Grab only the 2 columns we need
                pre_count_matrix = single_sample_counts[
                    ["Gene ID", "Counts Unstranded RNA-seq"]
                ]
                # Rename the counts column with the sample name
                pre_count_matrix = pre_count_matrix.rename(
                    columns={"Counts Unstranded RNA-seq": subdir[-11:]}
                )
                # Merge counts from each file
                if df.empty:
                    df = pre_count_matrix
                else:
                    df = df.merge(pre_count_matrix, how="outer", on="Gene ID")

    df = df.set_index("Gene ID").T

    return df


def sample_overview_data_wrangling(
    INPUT_PATH: str, DATA_TYPE: Literal["Bulk", "Dissociated Bulk"]
) -> pd.DataFrame:
    """From the samples_overview.xlsx spreadsheet extract the samples in each subset 
    and the site of origin, and clean up the dataframe.

    Output:
        Dataframe for either the bulk or dissociated bulk samples with columns for 
        site of origin and long sample id (includes which subset the data come from), 
        sorted by subset date and sample id.
    """

    if DATA_TYPE == "Bulk":
        cols = ["0509", "0626"]
    elif DATA_TYPE == "Dissociated Bulk":
        cols = ["0414", "0418"]

    df = pd.read_excel(INPUT_PATH + "samples_overview.xlsx", header=1, index_col=0)
    df = df[cols + ["Site of origin"]]
    df.drop(index=["Total across columns"], inplace=True)
    df = df.dropna(subset=cols, how="all")
    df = df.sort_index().sort_values(cols[0])

    index = []
    for col in cols:
        sub_df = df[df[col] == "X"]
        index += ["23" + col + "/" + str(x) for x in sub_df.index]

    df["Long Sample ID"] = index
    df = df[["Site of origin", "Long Sample ID"]]

    return df


def gene_ensemble_id_mapping(PATH: str) -> list[tuple[str, str]]:
    """Get mapping from ensemble gene ids to gene names.
    The GTF file was part of the reference genome download that we used for alignment.

    Output:
        list of tuples with the gene ID and gene name.
    """

    with open(PATH) as f:
        gtf = list(f)
        gtf = gtf[5:]

    gtf = [x for x in gtf if 'gene_id "' in x and 'gene_name "' in x]

    gtf = list(
        map(
            lambda x: (
                x.split('gene_id "')[1].split('"')[0],
                x.split('gene_name "')[1].split('"')[0],
            ),
            gtf,
        )
    )
    gtf = list(set(gtf))

    return gtf
