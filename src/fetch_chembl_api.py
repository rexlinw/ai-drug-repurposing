import requests
import pandas as pd
from tqdm import tqdm

BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"

def fetch_drugs(limit=100):
    url = f"{BASE_URL}/molecule.json?limit={limit}"
    response = requests.get(url)
    data = response.json()
    molecules = data['molecules']
    
    drug_data = []
    for mol in molecules:
        drug_data.append({
            "chembl_id": mol.get("molecule_chembl_id"),
            "name": mol.get("pref_name"),
            "max_phase": mol.get("max_phase"),
            "therapeutic_flag": mol.get("therapeutic_flag"),
            "molecule_type": mol.get("molecule_type")
        })
    
    return pd.DataFrame(drug_data)

def fetch_targets_for_drug(drug_id):
    url = f"{BASE_URL}/target.json?molecule_chembl_id={drug_id}"
    response = requests.get(url)
    data = response.json()
    return data.get('targets', [])

def fetch_activity_data(drug_id, limit=100):
    url = f"{BASE_URL}/activity.json?molecule_chembl_id={drug_id}&limit={limit}"
    response = requests.get(url)
    data = response.json()
    return data.get('activities', [])

if __name__ == "__main__":
    print("Fetching sample drug data from ChEMBL...")
    drugs_df = fetch_drugs(limit=100)
    drugs_df.to_csv("data/chembl_drugs_sample.csv", index=False)
    print("Saved sample drug data to data/chembl_drugs_sample.csv")
