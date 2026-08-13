"""
Inspect the scVI-integrated h5ad — verify what fields we actually have
before writing any visualization. Do NOT skip this step.
"""
import anndata as ad
import pandas as pd

path = "/Users/apple/Downloads/研究生毕业论文/Cross_Species/scvi_out/Fibroblast_scvi/integrated.h5ad"
adata = ad.read_h5ad(path)

print("=== SHAPE ===")
print(f"  cells   : {adata.n_obs}")
print(f"  genes   : {adata.n_vars}")
print()

print("=== .obs columns (per-cell metadata) ===")
for c in adata.obs.columns:
    dtype = adata.obs[c].dtype
    n_uniq = adata.obs[c].nunique()
    print(f"  {c:35s}  dtype={dtype}  n_unique={n_uniq}")
print()

print("=== .obsm keys (embeddings) ===")
for k in adata.obsm.keys():
    print(f"  {k:30s}  shape={adata.obsm[k].shape}")
print()

print("=== species distribution ===")
print(adata.obs["species"].value_counts())
print()

print("=== leiden distribution ===")
if "leiden_scvi" in adata.obs.columns:
    print(adata.obs["leiden_scvi"].value_counts().sort_index())
print()

# Cross-tab: species x leiden
print("=== species x leiden crosstab (helps decide whether integration mixed species) ===")
print(pd.crosstab(adata.obs["species"], adata.obs["leiden_scvi"]))
