import requests
import csv

def fetch_chembl_drugs(output_file="data/drugs.csv", limit=1000):
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule?limit={limit}&format=json"
    response = requests.get(url)

    if response.status_code == 200:
        drugs = response.json().get("molecules", [])
        with open(output_file, "w", newline="", encoding='utf-8') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["drug_id", "drug_name", "drug_type", "approval_status"])
            for d in drugs:
                phase = d.get("max_phase")
                if phase == 4:
                    status = "Approved"
                elif phase is None:
                    status = "Unknown"
                else:
                    status = "Investigational"

                writer.writerow([
                    d.get("molecule_chembl_id"),
                    d.get("pref_name") or "",
                    d.get("molecule_type") or "",
                    status
                ])
        print(f"Drug list saved to {output_file}")
    else:
        print("Failed to fetch drugs:", response.status_code)

fetch_chembl_drugs()
