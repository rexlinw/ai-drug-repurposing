import pandas as pd

def match_by_gene(drug_file, disease_file, output_file):
    """
    Match drugs to diseases based on shared gene targets.
    
    Args:
        drug_file (str): Path to the drugs.csv file.
        disease_file (str): Path to the diseases.csv file.
        output_file (str): Path to save the matched candidates.
    """
    print("[INFO] Loading input files...")
    drugs = pd.read_csv(drug_file)
    diseases = pd.read_csv(disease_file)

    print("[INFO] Matching drugs to diseases by gene...")
    matches = pd.merge(drugs, diseases, left_on="target_gene", right_on="gene_symbol")

    print(f"[INFO] Found {len(matches)} matched candidates.")
    matches.to_csv(output_file, index=False)
    print(f"[INFO] Saved matches to {output_file}")
