import pandas as pd

def merge_and_normalize(disease_file="data/diseases.csv", drug_file="data/drugs.csv", output_file="data/merged_candidates.csv"):
    diseases = pd.read_csv(disease_file)
    drugs = pd.read_csv(drug_file)

    diseases["disease_name"] = diseases["disease_name"].str.lower().str.strip()
    drugs["drug_name"] = drugs["drug_name"].str.lower().str.strip()

    merged = pd.merge(diseases.assign(key=1), drugs.assign(key=1), on="key").drop("key", axis=1)
    merged.to_csv(output_file, index=False)
    print(f"Merged candidates saved to {output_file}")

merge_and_normalize()
