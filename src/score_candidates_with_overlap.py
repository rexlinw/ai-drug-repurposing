import os
import pandas as pd

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(PROJECT_ROOT, "data")

drug_targets_df = pd.read_csv(os.path.join(DATA_DIR, "drug_gene_targets.csv"))
disease_signature_df = pd.read_csv(os.path.join(DATA_DIR, "disease_signature.csv"))

drug_targets_df = drug_targets_df.dropna(subset=["gene"])
drug_targets_df["gene"] = drug_targets_df["gene"].astype(str).str.upper()
disease_genes = set(disease_signature_df["gene_symbol"].dropna().astype(str).str.upper())

drug_to_genes = drug_targets_df.groupby("chembl_id")["gene"].apply(set).to_dict()
scored = []
for chembl_id, gene_set in drug_to_genes.items():
    overlap = gene_set & disease_genes
    score = len(overlap)
    scored.append({
        "chembl_id": chembl_id,
        "overlap_score": score,
        "num_drug_genes": len(gene_set),
        "num_overlap_genes": len(overlap),
        "overlap_genes": ", ".join(overlap)
    })

scored_df = pd.DataFrame(scored).sort_values(by="overlap_score", ascending=False)
scored_df.to_csv(os.path.join(DATA_DIR, "scored_candidates_overlap.csv"), index=False)
print("✅ Saved: scored_candidates_overlap.csv")
