import requests
import csv
import time

gene_list = ["5742", "2157", "5290"]
KEGG_BASE = "https://rest.kegg.jp"

def fetch_pathways(entrez_id):
    url = f"{KEGG_BASE}/link/pathway/hsa:{entrez_id}"
    response = requests.get(url)
    if response.status_code != 200:
        print(f"Failed to fetch pathways for gene {entrez_id}")
        return []
    lines = response.text.strip().split("\n")
    pathways = []
    for line in lines:
        parts = line.split("\t")
        if len(parts) == 2:
            pathways.append(parts[1])
    return pathways

def main():
    with open("data/kegg_gene_pathways.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["entrez_id", "pathway_id"])

        for gene in gene_list:
            print(f"Fetching pathways for gene {gene}")
            pathways = fetch_pathways(gene)
            for p in pathways:
                writer.writerow([gene, p])
            time.sleep(1)

    print("Saved kegg_gene_pathways.csv")

if __name__ == "__main__":
    main()
