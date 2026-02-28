import os
import pandas as pd

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(PROJECT_ROOT, "data")

drug_df = pd.read_csv(os.path.join(DATA_DIR, "drug_gene_targets.csv"))
print("Drug target columns:", drug_df.columns)
disease_df = pd.read_csv(os.path.join(DATA_DIR, "disease_signature.csv"))
disease_df["entrez_id"] = disease_df["entrez_id"].astype(str)

kegg_df = pd.read_csv(os.path.join(DATA_DIR, "kegg_gene_pathways.csv"))
kegg_df["entrez_id"] = kegg_df["entrez_id"].astype(str)

gene_to_pathways = kegg_df.groupby("entrez_id")["pathway_id"].apply(set).to_dict()
disease_pathways = set()
for gene in disease_df["entrez_id"]:
    disease_pathways |= gene_to_pathways.get(str(gene), set())

drug_scores = []
for chembl_id, group in drug_df.groupby("chembl_id"):
    drug_genes = group["gene"]
    drug_pathways = set()
    for gene in drug_genes:
        drug_pathways |= gene_to_pathways.get(str(gene), set())
    
    overlap = len(drug_pathways & disease_pathways)
    union = len(drug_pathways | disease_pathways)
    jaccard = overlap / union if union != 0 else 0

    drug_scores.append({
        "chembl_id": chembl_id,
        "pathway_overlap_score": jaccard
    })

pd.DataFrame(drug_scores).to_csv(os.path.join(DATA_DIR, "scored_candidates_pathway.csv"), index=False)
print("✅ Done: scored_candidates_pathway.csv saved.")

