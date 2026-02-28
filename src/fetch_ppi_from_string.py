import requests
import pandas as pd
import os

def fetch_string_ppi(genes, species=9606, score_threshold=900):
    """
    Fetch high-confidence PPI data from STRING for a list of human genes.
    """
    print(f"Fetching PPI data from STRING for {len(genes)} genes...")
    
    url = "https://string-db.org/api/tsv/network"
    identifiers = "%0d".join(genes)
    params = {
        "identifiers": identifiers,
        "species": species,
        "required_score": score_threshold,
        "limit": 1000,
        "caller_identity": "ai-drug-repurposing"
    }

    response = requests.post(url, data=params)

    if response.status_code != 200:
        raise Exception(f"Error from STRING API: {response.status_code}\n{response.text}")

    lines = response.text.strip().split("\n")
    header = lines[0].split("\t")
    data = [line.split("\t") for line in lines[1:]]
    df = pd.DataFrame(data, columns=header)

    return df


def main():
    gene_list = [
        "TP53", "EGFR", "BRCA1", "TNF", "MTOR", "IL6", "VEGFA", "AKT1",
        "JUN", "MYC", "MMP9", "INS", "CXCL8", "CDK2"
    ]
    df = fetch_string_ppi(gene_list)
    os.makedirs("data/ppi", exist_ok=True)
    output_path = "data/ppi/string_ppi_subset.csv"
    df.to_csv(output_path, index=False)
    print(f"\n✅ Saved STRING PPI subset to {output_path}")

if __name__ == "__main__":
    main()
