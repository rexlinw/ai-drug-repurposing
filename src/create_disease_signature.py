import pandas as pd

genes = ['TP53', 'EGFR', 'TNF', 'IL6', 'VEGFA']
df = pd.DataFrame({'gene_symbol': genes})
df.to_csv('../data/disease_signature.csv', index=False)

print("✅ disease_signature.csv created in ../data/")
