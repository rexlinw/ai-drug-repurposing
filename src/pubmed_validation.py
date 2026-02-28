import os
from Bio import Entrez
import time

Entrez.email = os.environ.get("NCBI_EMAIL", "your-email@example.com")

def pubmed_hit_count(drug_name, disease_name, delay=0.3):
    """
    Searches PubMed for articles that mention both drug and disease.
    
    Args:
        drug_name (str): Name of the drug.
        disease_name (str): Name of the disease.
        delay (float): Delay between requests to avoid rate limits.

    Returns:
        int: Number of matching PubMed articles.
    """
    query = f'"{drug_name}" AND "{disease_name}"'
    try:
        handle = Entrez.esearch(db="pubmed", term=query)
        record = Entrez.read(handle)
        time.sleep(delay)
        return int(record["Count"])
    except Exception as e:
        print(f"[ERROR] PubMed query failed for {query}: {e}")
        return 0
