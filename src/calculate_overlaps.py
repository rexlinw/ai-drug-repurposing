import csv
from collections import defaultdict

def load_disease_genes(filename):
    disease_genes = {}
    with open(filename) as f:
        reader = csv.DictReader(f)
        for row in reader:
            genes = row["entrez_genes"].split(",") if row["entrez_genes"] else []
            disease_genes[row["disease_name"]] = set(genes)
    return disease_genes

def load_drug_genes(filename):
    drug_genes = {}
    with open(filename) as f:
        reader = csv.DictReader(f)
        for row in reader:
            genes = row["entrez_genes"].split(",") if row["entrez_genes"] else []
            drug_genes[row["drug_name"]] = set(genes)
    return drug_genes

def load_gene_pathways(filename):
    gene_pathways = defaultdict(set)
    with open(filename) as f:
        reader = csv.DictReader(f)
        for row in reader:
            gene_pathways[row["entrez_id"]].add(row["pathway_id"])
    return gene_pathways

def pathway_overlap_score(genes1, genes2, gene_pathways):
    pathways1 = set()
    pathways2 = set()
    for g in genes1:
        pathways1.update(gene_pathways.get(g, []))
    for g in genes2:
        pathways2.update(gene_pathways.get(g, []))
    intersection = pathways1.intersection(pathways2)
    union = pathways1.union(pathways2)
    if not union:
        return 0.0
    return len(intersection) / len(union)

def main():
    disease_genes = load_disease_genes("data/disease_gene_map.csv")
    drug_genes = load_drug_genes("data/drug_target_genes.csv")
    gene_pathways = load_gene_pathways("data/kegg_gene_pathways.csv")

    results = []
    for disease, dgenes in disease_genes.items():
        for drug, dggenes in drug_genes.items():
            score = pathway_overlap_score(dgenes, dggenes, gene_pathways)
            if score > 0:
                results.append((disease, drug, score))

    results.sort(key=lambda x: x[2], reverse=True)
    with open("data/disease_drug_pathway_overlap.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["disease", "drug", "pathway_overlap_score"])
        for row in results:
            writer.writerow(row)

    print(f"Done! Results saved in data/disease_drug_pathway_overlap.csv")

if __name__ == "__main__":
    main()
