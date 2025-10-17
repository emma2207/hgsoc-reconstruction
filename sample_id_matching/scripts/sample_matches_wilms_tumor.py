import pandas as pd

if __name__ == "__main__":
    df = pd.read_excel("../data/sample_matches/SraRunTable_sn_wilms_tumor_AJMupdated.xlsx")

    df = df[["Sample Name", "Bulk RNA-seq sample ID"]]
    df.rename(columns={
        "Sample Name": "Single-nucleus sample ID", 
        "Bulk RNA-seq sample ID": "Bulk sample ID"
        }, inplace=True)

    print(df.head())
    df.to_csv("../data/sample_matches/sample_matches_wilms_tumor.csv", index=False)