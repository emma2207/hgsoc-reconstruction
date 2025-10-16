import pandas as pd

if __name__ == "__main__":
    df = pd.read_excel("../data/sample_matches/sample_matches_central_nervous_system_tumor.xlsx")

    temp_df = df["Pool and multiplex ID"].str.split("-", n=1, expand=True)
    temp_df.columns = ["Pool", "Multiplex ID"]
    df = pd.concat([df, temp_df], axis=1)
    df.drop(columns=["Pool and multiplex ID"], inplace=True)

    # Isolate and clean single cell samples
    single_cell_df = df[df["Pool"].str.contains("Sample")]
    single_cell_df = single_cell_df.dropna(axis=1, how="all")
    single_cell_df.rename(columns={"Pool": "Sample ID"}, inplace=True)
    single_cell_df["Sample ID"] = single_cell_df["Sample ID"].str.strip()
    # single_cell_df.to_csv("../data/sample_matches/sample_matches_central_nervous_system_tumor.csv", index=False)

    # Isolate and clean bulk samples
    bulk_df = df[df["Pool"].str.contains("Pool")]
    # Split HTO derived cDNA derived cDNA into 3 columns
    temp_df = bulk_df["HTO derived cDNA"].str.split("/", n=2, expand=True)
    temp_df.columns = ["HTO derived cDNA - 1", "HTO derived cDNA - 2", "HTO derived cDNA - 3"]
    bulk_df = pd.concat([bulk_df, temp_df], axis=1)
    bulk_df.drop(columns=["HTO derived cDNA"], inplace=True)
    # Split mRNA derived cDNA derived cDNA into 3 columns
    temp_df = bulk_df["mRNA derived cDNA"].str.split("/", n=2, expand=True)
    temp_df.columns = ["mRNA derived cDNA - 1", "mRNA derived cDNA - 2", "mRNA derived cDNA - 3"]
    bulk_df = pd.concat([bulk_df, temp_df], axis=1)
    bulk_df.drop(columns=["mRNA derived cDNA"], inplace=True)
    # Split Pool into 3 columns
    temp_df = bulk_df["Pool"].str.split(", ", n=2, expand=True)
    temp_df.columns = ["Pool - 1", "Pool - 2", "Pool - 3"]
    bulk_df = pd.concat([bulk_df, temp_df], axis=1)
    bulk_df.drop(columns=["Pool"], inplace=True)

    # Unpivot wide dataframe
    long_df = pd.wide_to_long(bulk_df, stubnames=["HTO derived cDNA - ", "mRNA derived cDNA - ", "Pool - "], i="Name", j="Sample Number", sep="")
    long_df.columns = long_df.columns.str.rstrip(" - ")
    long_df = long_df.dropna(subset=["HTO derived cDNA", "mRNA derived cDNA", "Pool"], how="all")
    long_df["Pool"] = long_df["Pool"].str.strip()
    long_df = long_df.sort_values(by=["Name", "Sample Number"])

    long_df.to_csv("../data/sample_matches/sample_matches_cnst_long.csv", index=False)

    print(long_df.head(10))