import pandas as pd

pl_path = "/pl/active/cgreene-sc-hgsoc/mismatch_project_data/metadata"

if __name__ == "__main__":
    df = pd.read_excel(
        pl_path
        + "/received_match_info/sample_matches_central_nervous_system_tumor.xlsx"
    )

    # Split Pool and multiplex ID into 2 columns
    temp_df = df["Pool and multiplex ID"].str.split("-", n=1, expand=True)
    temp_df.columns = ["Pool", "Multiplex ID"]
    df = pd.concat([df, temp_df], axis=1)
    df.drop(columns=["Pool and multiplex ID"], inplace=True)

    # Isolate and clean non-pooled samples
    non_pooled_df = df[df["Pool"].str.contains("Sample")]
    non_pooled_df = non_pooled_df.dropna(axis=1, how="all")
    non_pooled_df.rename(columns={"Pool": "Sample ID"}, inplace=True)
    non_pooled_df["Sample ID"] = non_pooled_df["Sample ID"].str.strip()

    # Isolate and clean pooled samples
    pooled_df = df[df["Pool"].str.contains("Pool")]
    # Split HTO derived cDNA derived cDNA into 3 columns
    temp_df = pooled_df["HTO derived cDNA"].str.split("/", n=2, expand=True)
    temp_df.columns = [
        "HTO derived cDNA - 1",
        "HTO derived cDNA - 2",
        "HTO derived cDNA - 3",
    ]
    pooled_df = pd.concat([pooled_df, temp_df], axis=1)
    pooled_df.drop(columns=["HTO derived cDNA"], inplace=True)
    # Split mRNA derived cDNA derived cDNA into 3 columns
    temp_df = pooled_df["mRNA derived cDNA"].str.split("/", n=2, expand=True)
    temp_df.columns = [
        "mRNA derived cDNA - 1",
        "mRNA derived cDNA - 2",
        "mRNA derived cDNA - 3",
    ]
    pooled_df = pd.concat([pooled_df, temp_df], axis=1)
    pooled_df.drop(columns=["mRNA derived cDNA"], inplace=True)
    # Split Pool into 3 columns
    temp_df = pooled_df["Pool"].str.split(", ", n=2, expand=True)
    temp_df.columns = ["Pool - 1", "Pool - 2", "Pool - 3"]
    pooled_df = pd.concat([pooled_df, temp_df], axis=1)
    pooled_df.drop(columns=["Pool"], inplace=True)

    # Unpivot wide dataframe
    long_df = pd.wide_to_long(
        pooled_df,
        stubnames=["HTO derived cDNA - ", "mRNA derived cDNA - ", "Pool - "],
        i="Name",
        j="Sample Number",
        sep="",
    )
    long_df.columns = long_df.columns.str.rstrip(" - ")
    long_df = long_df.dropna(
        subset=["HTO derived cDNA", "mRNA derived cDNA", "Pool"], how="all"
    )
    long_df["Pool"] = long_df["Pool"].str.strip()
    long_df = long_df.sort_values(by=["Name", "Sample Number"])

    # Get pool contents and unpivot
    pool_contents_df = long_df[["Pool", "Bulk RNA seq"]].copy()
    pool_contents_df["Pool"] = (
        pool_contents_df["Pool"].str.replace("Pool", "").str.replace("_", "")
    )
    pool_contents_df = (
        pool_contents_df.drop_duplicates(ignore_index=True)
        .sort_values(by=["Pool"])
        .reset_index(drop=True)
    )
    pool_contents_df["idx"] = pool_contents_df.groupby("Pool").cumcount()
    pivoted_pool_contents_df = pool_contents_df.pivot(
        index="Pool", columns="idx", values="Bulk RNA seq"
    )
    pivoted_pool_contents_df = pivoted_pool_contents_df.dropna(how="all", axis=1)

    pivoted_pool_contents_df.to_csv(
        pl_path + "/pool_contents/pool_contents_central_nervous_system_tumor.csv",
        index=True,
    )

    # Combine sample matches
    long_df = long_df[["Bulk RNA seq", "HTO derived cDNA", "mRNA derived cDNA"]]
    long_df.index = pd.Index([idx[0] for idx in list(long_df.index)])
    long_df = long_df.drop_duplicates(ignore_index=False)
    mapping_df = pd.concat(
        [
            long_df,
            non_pooled_df[["Name", "Bulk RNA seq", "mRNA derived cDNA"]].set_index(
                "Name"
            ),
        ],
        axis=0,
    )

    mapping_df.to_csv(
        pl_path + "/sample_matches/sample_matches_central_nervous_system_tumor.csv",
        index=True,
    )
