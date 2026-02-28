import argparse
import os
import requests
import csv
import sys
import time
import pandas as pd

PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

DATA_DIR = os.path.join(PROJECT_ROOT, "data")
os.makedirs(DATA_DIR, exist_ok=True)

def save_dicts_to_csv(data, keys, filename):
    with open(filename, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        for row in data:
            writer.writerow({k: row.get(k, "") for k in keys})
    print(f"Saved {filename}")

def fetch_disgenet_diseases():
    print("Fetching diseases from DisGeNET...")
    api_key = os.environ.get("DISGENET_API_KEY", "")

    if api_key:
        headers = {"Accept": "application/json", "Authorization": f"Bearer {api_key}"}
        response = requests.get("https://www.disgenet.org/api/v1/diseases", headers=headers, timeout=30)
        if response.status_code == 200:
            try:
                data = response.json()
                diseases = data if isinstance(data, list) else data.get("results", data.get("data", []))
                if diseases:
                    keys = ["diseaseid", "disease_name", "disease_type"] if diseases else []
                    if diseases and isinstance(diseases[0], dict):
                        all_keys = list(diseases[0].keys())
                        mapping = {"diseaseId": "diseaseid", "name": "disease_name", "type": "disease_type"}
                        rows = []
                        for d in diseases[:500]:
                            rows.append({
                                "diseaseid": d.get("diseaseId") or d.get("diseaseid") or d.get("umlsId", ""),
                                "disease_name": d.get("name") or d.get("disease_name", ""),
                                "disease_type": d.get("type") or d.get("disease_type", "disease")
                            })
                        save_dicts_to_csv(rows, ["diseaseid", "disease_name", "disease_type"], os.path.join(DATA_DIR, "diseases.csv"))
                    else:
                        save_dicts_to_csv(diseases, keys, os.path.join(DATA_DIR, "diseases.csv"))
                    return
            except Exception as e:
                print(f"DisGeNET v1 parse error: {e}")

    print("Using seed disease data (set DISGENET_API_KEY for live API).")
    seed_diseases = [
        {"diseaseid": "C0002395", "disease_name": "Alzheimer's disease", "disease_type": "disease"},
        {"diseaseid": "C0006142", "disease_name": "Breast cancer", "disease_type": "disease"},
        {"diseaseid": "C0011849", "disease_name": "Diabetes mellitus", "disease_type": "disease"},
        {"diseaseid": "C0021368", "disease_name": "Inflammation", "disease_type": "disease"},
    ]
    save_dicts_to_csv(seed_diseases, ["diseaseid", "disease_name", "disease_type"], os.path.join(DATA_DIR, "diseases.csv"))

def fetch_chembl_drugs(limit=500):
    print("Fetching drugs from ChEMBL...")
    base_url = "https://www.ebi.ac.uk/chembl/api/data/molecule"
    all_drugs = []
    page_size = 100
    for offset in range(0, limit, page_size):
        url = f"{base_url}?limit={page_size}&offset={offset}&format=json&max_phase=4"
        resp = requests.get(url)
        if resp.status_code != 200:
            print(f"Warning: failed at offset {offset}")
            break
        data = resp.json()
        molecules = data.get("molecules", [])
        if not molecules:
            break
        all_drugs.extend(molecules)
        time.sleep(0.2)
    drugs_clean = []
    for d in all_drugs:
        phase = d.get("max_phase")
        status = "Approved" if phase == 4 else ("Investigational" if phase else "Unknown")
        drugs_clean.append({
            "drug_id": d.get("molecule_chembl_id"),
            "drug_name": d.get("pref_name") or "",
            "drug_type": d.get("molecule_type") or "",
            "approval_status": status
        })
    keys = ["drug_id", "drug_name", "drug_type", "approval_status"]
    save_dicts_to_csv(drugs_clean, keys, os.path.join(DATA_DIR, "drugs.csv"))

def merge_candidates():
    print("Merging disease and drug candidates...")
    diseases = pd.read_csv(os.path.join(DATA_DIR, "diseases.csv"))
    drugs = pd.read_csv(os.path.join(DATA_DIR, "drugs.csv"))
    diseases["disease_name"] = diseases["disease_name"].str.lower().str.strip()
    drugs["drug_name"] = drugs["drug_name"].str.lower().str.strip()

    diseases["key"] = 1
    drugs["key"] = 1
    merged = pd.merge(diseases, drugs, on="key").drop("key", axis=1)
    merged.to_csv(os.path.join(DATA_DIR, "merged_candidates.csv"), index=False)
    print(f"Merged candidates saved to {os.path.join(DATA_DIR, 'merged_candidates.csv')}")

def fetch_gene_disease(limit=500):
    print("Fetching gene-disease associations from DisGeNET...")
    api_key = os.environ.get("DISGENET_API_KEY", "")
    headers = {"Accept": "application/json"}
    if api_key:
        headers["Authorization"] = f"Bearer {api_key}"

    response = requests.get(
        "https://www.disgenet.org/api/gda/disease",
        headers=headers,
        params={"limit": limit},
        timeout=30
    )
    if response.status_code == 200:
        try:
            gda = response.json()
            records = [{"disease_id": e.get("diseaseid"), "gene_id": str(e.get("geneid", "")), "score": e.get("score")} for e in gda]
            if records:
                save_dicts_to_csv(records, ["disease_id", "gene_id", "score"], os.path.join(DATA_DIR, "gene_disease.csv"))
                return
        except Exception as e:
            print(f"Gene-disease parse error: {e}")

    print("Using seed gene-disease data.")
    disease_genes = {
        "C0002395": ["APOE", "PSEN1", "APP", "TNF", "IL6"],
        "C0006142": ["BRCA1", "BRCA2", "TP53", "EGFR", "VEGFA"],
        "C0011849": ["INS", "INSR", "TNF", "IL6", "VEGFA"],
        "C0021368": ["PTGS1", "PTGS2", "TNF", "IL6"],
    }
    records = []
    for did, genes in disease_genes.items():
        for g in genes:
            records.append({"disease_id": did, "gene_id": g, "score": 0.5})
    save_dicts_to_csv(records, ["disease_id", "gene_id", "score"], os.path.join(DATA_DIR, "gene_disease.csv"))

def fetch_gene_drug(limit=100):
    drug_targets_path = os.path.join(DATA_DIR, "drug_gene_targets.csv")
    if os.path.exists(drug_targets_path):
        print("Using existing drug_gene_targets.csv.")
        return
    from src.fetch_drug_gene_targets import run_fetch_drug_gene_targets
    run_fetch_drug_gene_targets(limit=limit, data_dir=DATA_DIR)

def fetch_ppi_pathways():
    print("PPI and pathway data are large. Download from STRING, KEGG manually.")
    print("Then parse and save as CSVs for gene interactions and pathways.")

def calculate_overlap_scores():
    print("Calculating overlap scores...")
    try:
        gene_disease = pd.read_csv(os.path.join(DATA_DIR, "gene_disease.csv"))
        merged = pd.read_csv(os.path.join(DATA_DIR, "merged_candidates.csv"))

        drug_targets_path = os.path.join(DATA_DIR, "drug_gene_targets.csv")
        known_drug_genes = {
            "CHEMBL25": ["PTGS1", "PTGS2"],
            "CHEMBL1431": ["INSR", "AMPK"],
            "CHEMBL198": ["EGFR"],
            "CHEMBL1201582": ["TNF"],
            "CHEMBL472": ["VEGFA"],
            "CHEMBL1618": ["IL6"],
        }
        if os.path.exists(drug_targets_path):
            dt = pd.read_csv(drug_targets_path)
            if "chembl_id" in dt.columns and "gene" in dt.columns:
                dt["gene"] = dt["gene"].astype(str).str.upper()
                for cid, grp in dt.groupby("chembl_id"):
                    genes = set(grp["gene"].dropna())
                    known_drug_genes[str(cid)] = list(genes) if genes else known_drug_genes.get(str(cid), [])

        drug_genes = {k: set(v) for k, v in known_drug_genes.items()}
        disease_genes = gene_disease.groupby("disease_id")["gene_id"].apply(lambda x: set(str(g).upper() for g in x)).to_dict()

        scores = []
        for idx, row in merged.iterrows():
            d_id = row["diseaseid"] if "diseaseid" in row else row["disease_id"]
            dr_id = row["drug_id"]
            d_genes = disease_genes.get(d_id, set())
            dr_genes = drug_genes.get(dr_id, set())
            if d_genes and dr_genes:
                intersection = d_genes.intersection(dr_genes)
                union = d_genes.union(dr_genes)
                jaccard = len(intersection) / len(union) if union else 0
                scores.append([d_id, row["disease_name"], dr_id, row["drug_name"], jaccard])
            else:
                scores.append([d_id, row["disease_name"], dr_id, row["drug_name"], 0])
        scored_df = pd.DataFrame(scores, columns=["disease_id", "disease_name", "drug_id", "drug_name", "overlap_score"])
        scored_df.sort_values(by="overlap_score", ascending=False, inplace=True)
        scored_df.to_csv(os.path.join(DATA_DIR, "scored_candidates.csv"), index=False)
        print(f"Scored candidates saved to {os.path.join(DATA_DIR, 'scored_candidates.csv')}")
    except Exception as e:
        print(f"Error calculating scores: {e}")

def validate_pubmed(top_n=50):
    """Add PubMed co-mention counts to top N scored candidates. Requires NCBI_EMAIL env var."""
    from src.pubmed_validation import pubmed_hit_count

    email = os.environ.get("NCBI_EMAIL", "").strip()
    if not email or email == "your-email@example.com":
        print("Skipping PubMed validation: set NCBI_EMAIL env var (e.g. export NCBI_EMAIL=you@example.com)")
        return

    scored_path = os.path.join(DATA_DIR, "scored_candidates.csv")
    if not os.path.exists(scored_path):
        print("No scored_candidates.csv found. Run pipeline without --pubmed first.")
        return

    print(f"Validating top {top_n} candidates with PubMed...")
    df = pd.read_csv(scored_path).head(top_n)
    counts = []
    for idx, row in df.iterrows():
        drug_name = str(row.get("drug_name", "")).strip()
        disease_name = str(row.get("disease_name", "")).strip()
        if drug_name and disease_name:
            c = pubmed_hit_count(drug_name, disease_name, delay=0.35)
            counts.append(c)
        else:
            counts.append(0)
    df["pubmed_count"] = counts
    out_path = os.path.join(DATA_DIR, "scored_candidates_validated.csv")
    df.to_csv(out_path, index=False)
    print(f"Saved {out_path}")

def main():
    parser = argparse.ArgumentParser(description="AI Drug Repurposing Pipeline")
    parser.add_argument("--pubmed", action="store_true", help="Add PubMed validation to top candidates")
    parser.add_argument("--top-n", type=int, default=50, help="Number of top candidates to validate with PubMed (default: 50)")
    args = parser.parse_args()

    fetch_disgenet_diseases()
    fetch_chembl_drugs(limit=500)
    merge_candidates()
    fetch_gene_disease(limit=500)
    fetch_gene_drug(limit=100)
    fetch_ppi_pathways()
    calculate_overlap_scores()

    if args.pubmed:
        validate_pubmed(top_n=args.top_n)

if __name__ == "__main__":
    main()
