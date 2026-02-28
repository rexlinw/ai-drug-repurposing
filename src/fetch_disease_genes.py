import requests
import csv
import time

diseases = {
    "Alzheimer's disease": "C0002395",
    "Breast cancer": "C0006142",
    "Diabetes mellitus": "C0011849",
}

API_BASE = "https://www.disgenet.org/api/gda/disease/"

def fetch_genes_for_disease(disease_id):
    url = f"{API_BASE}{disease_id}?source=ALL"
    headers = {"Accept": "application/json"}
    response = requests.get(url, headers=headers)

    if response.status_code != 200:
        print(f"Failed to fetch {disease_id}: {response.status_code}")
        return []

    data = response.json()
    genes = set()
    for record in data:
        geneid = record.get("geneid")
        if geneid:
            genes.add(str(geneid))
    return list(genes)

def main():
    with open("data/disease_gene_map.csv", "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["disease_name", "disease_id", "entrez_genes"])

        for name, d_id in diseases.items():
            print(f"Fetching genes for {name} ({d_id})...")
            genes = fetch_genes_for_disease(d_id)
            if genes:
                writer.writerow([name, d_id, ",".join(genes)])
            else:
                writer.writerow([name, d_id, ""])
            time.sleep(1)

    print("Done! disease_gene_map.csv saved in data/")

if __name__ == "__main__":
    main()
