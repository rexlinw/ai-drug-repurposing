import os
import pandas as pd
import requests
import time

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(PROJECT_ROOT, "data")

def get_targets_for_drug(chembl_id):
    url = f"https://www.ebi.ac.uk/chembl/api/data/target?molecule_chembl_id={chembl_id}&format=json"
    try:
        response = requests.get(url)
        if response.status_code == 200:
            return response.json()["targets"]
        else:
            return []
    except:
        return []

def get_gene_symbol_for_target(target_id):
    url = f"https://www.ebi.ac.uk/chembl/api/data/target/{target_id}.json"
    try:
        response = requests.get(url)
        if response.status_code == 200:
            target = response.json()
            for xref in target.get("target_components", [{}])[0].get("target_component_xrefs", []):
                if "xref_id" in xref:
                    return xref["xref_id"]
        return None
    except:
        return None

def run_fetch_drug_gene_targets(limit=100, data_dir=None):
    """Fetch drug-gene targets from ChEMBL. Uses drugs.csv or chembl_drugs_sample.csv."""
    data_dir = data_dir or DATA_DIR
    drugs_csv = os.path.join(data_dir, "drugs.csv")
    chembl_sample_csv = os.path.join(data_dir, "chembl_drugs_sample.csv")

    if os.path.exists(drugs_csv):
        drug_df = pd.read_csv(drugs_csv)
        chembl_col = "drug_id" if "drug_id" in drug_df.columns else "chembl_id"
    elif os.path.exists(chembl_sample_csv):
        drug_df = pd.read_csv(chembl_sample_csv)
        chembl_col = "chembl_id"
    else:
        print("No drugs.csv or chembl_drugs_sample.csv found. Skipping drug-gene fetch.")
        return

    chembl_ids = drug_df[chembl_col].dropna().unique().tolist()[:limit]
    print(f"Fetching targets for {len(chembl_ids)} drugs from ChEMBL...")

    results = []
    for idx, chembl_id in enumerate(chembl_ids):
        targets = get_targets_for_drug(chembl_id)
        for target in targets:
            target_id = target.get("target_chembl_id")
            if target_id:
                gene = get_gene_symbol_for_target(target_id)
                if gene:
                    results.append({"chembl_id": chembl_id, "target_id": target_id, "gene": gene})
        time.sleep(0.5)
        if (idx + 1) % 20 == 0:
            print(f"  Processed {idx + 1}/{len(chembl_ids)}...")

    if results:
        df = pd.DataFrame(results)
        df.to_csv(os.path.join(data_dir, "drug_gene_targets.csv"), index=False)
        print(f"Saved drug_gene_targets.csv ({len(results)} mappings)")
    else:
        print("No drug-gene mappings collected.")

if __name__ == "__main__":
    run_fetch_drug_gene_targets(limit=100)



