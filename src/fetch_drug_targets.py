import requests
import csv
import time

drug_list = {
    "Aspirin": "CHEMBL25",
    "Metformin": "CHEMBL1431",
    "Atorvastatin": "CHEMBL1487",
}

def fetch_targets(chembl_id):
    url = f"https://www.ebi.ac.uk/chembl/api/data/target?molecule_chembl_id={chembl_id}&format=json"
    response = requests.get(url)
    if response.status_code != 200:
        print(f"Failed to fetch targets for {chembl_id}")
        return []
    data = response.json()
    targets = set()
    for target in data.get("targets", []):
        if "target_components" in target:
            for comp in target["target_components"]:
                gene_id = comp.get("gene_id")
                if gene_id:
                    targets.add(str(gene_id))
    return list(targets)

def main():
    with open("data/drug_target_genes.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["drug_name", "chembl_id", "entrez_genes"])

        for name, chembl_id in drug_list.items():
            print(f"Fetching targets for {name} ({chembl_id})")
            genes = fetch_targets(chembl_id)
            writer.writerow([name, chembl_id, ",".join(genes)])
            time.sleep(1)

    print("Saved drug_target_genes.csv")

if __name__ == "__main__":
    main()
