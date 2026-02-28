## AI-Driven Drug Repurposing Pipeline

This repository implements a small, end‑to‑end **drug repurposing pipeline** that combines:
- Public biomedical APIs (DisGeNET, ChEMBL)
- Simple gene / pathway overlap scoring
- Optional PubMed literature validation

The goal is to generate candidate **disease–drug pairs** where the drug’s target genes overlap with genes implicated in a disease, then optionally check whether the pair is discussed in the literature.

---

## Project structure

- `run_pipeline.py` – main entry point; orchestrates data fetching and overlap scoring
- `src/fetch_drug_gene_targets.py` – fetch drug → target gene mappings from ChEMBL
- `src/score_candidates_with_overlap.py` – gene‑level overlap scoring utility
- `src/score_pathway_overlap.py` – pathway‑level overlap scoring utility (requires KEGG + mapped IDs)
- `src/pubmed_validation.py` – PubMed co‑mention counting for candidate pairs
- `data/` – all intermediate and output CSV files

---

## Setup

From the project root:

```bash
python3 -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

---

## Environment variables

These are optional but recommended:

- **`DISGENET_API_KEY`** – API token from DisGeNET (free registration).  
  If unset, the pipeline falls back to a small set of seed diseases and gene–disease associations.

- **`NCBI_EMAIL`** – your email address for NCBI Entrez (PubMed) usage.  
  Required if you enable PubMed validation so that requests comply with NCBI’s usage policy.

Example:

```bash
export DISGENET_API_KEY="your_disgenet_token"
export NCBI_EMAIL="you@example.com"
```

---

## Running the core pipeline

Basic run (fetch diseases/drugs, merge candidates, compute gene‑overlap scores):

```bash
source .venv/bin/activate
python run_pipeline.py
```

Key outputs in `data/`:

- `diseases.csv` – diseases (from DisGeNET or seed set)
- `drugs.csv` – approved drugs from ChEMBL
- `merged_candidates.csv` – Cartesian product of diseases × drugs
- `gene_disease.csv` – gene–disease associations (from DisGeNET or seed set)
- `scored_candidates.csv` – disease–drug pairs ranked by gene‑set overlap

---

## Optional: PubMed validation

When enabled, the pipeline can re‑score the **top N candidates** using PubMed co‑mentions of `(drug_name AND disease_name)`, and write:

- `data/scored_candidates_validated.csv` – same rows as `scored_candidates.csv` plus a `pubmed_count` column.

Usage:

```bash
export NCBI_EMAIL="you@example.com"
python run_pipeline.py --pubmed --top-n 50
```

---

## Optional utilities

- `src/fetch_drug_gene_targets.py`  
  Uses ChEMBL to build `data/drug_gene_targets.csv` with `chembl_id → gene` mappings.  
  By default it looks for `data/chembl_drugs_sample.csv`; we will also align it to use `data/drugs.csv` to stay consistent with the main pipeline.

- `src/score_candidates_with_overlap.py`  
  Standalone script to score drugs against a disease gene signature in `data/disease_signature.csv` and write an additional scored CSV.

- `src/score_pathway_overlap.py`  
  Demonstrates pathway‑level scoring using `data/kegg_gene_pathways.csv`. It assumes consistent Entrez IDs in both the disease signature and KEGG file; if your signature is gene symbols, you will need an ID‑mapping step first.

---

## Next steps / extensions

- Improve pathway scoring by adding robust gene‑ID mapping (symbol ↔ Entrez)
- Add ML‑based models that combine chemical descriptors (RDKit), overlap scores, and PubMed features
- Build a small API or UI to explore ranked candidate pairs interactively

