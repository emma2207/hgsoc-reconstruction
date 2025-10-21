import pandas as pd

pl_path = "pl/active/cgreene-sc-hgsoc/mismatch_project_data/metadata"

if __name__ == "__main__":
    df = pd.read_excel(
        pl_path + "/received_match_info/SraRunTable_sn_wilms_tumor_AJMupdated.xlsx"
    )

    df = df[["Sample Name", "Bulk RNA-seq sample ID"]]
    df.rename(
        columns={
            "Sample Name": "Single-nucleus sample ID",
            "Bulk RNA-seq sample ID": "Bulk sample ID",
        },
        inplace=True,
    )

    df.to_csv(pl_path + "/sample_matches/sample_matches_wilms_tumor.csv", index=False)
