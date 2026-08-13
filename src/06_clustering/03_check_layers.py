"""
05a — Before running marker analysis, verify what's actually stored
in the integrated.h5ad:
  - .X            main matrix (raw counts? log-norm? scVI latent?)
  - .raw          typically holds pre-integration normalized data
  - .layers       any additional matrices (counts, log1p, etc.)
  - .obsm         embeddings (X_scVI, X_umap)
"""
import anndata as ad
import numpy as np

path = "/Users/apple/Downloads/研究生毕业论文/Cross_Species/scvi_out/Fibroblast_scvi/integrated.h5ad"
adata = ad.read_h5ad(path)

print("=== SHAPE ===")
print(f"  {adata.n_obs} cells × {adata.n_vars} genes\n")

print("=== .X summary ===")
X = adata.X
print(f"  dtype  : {X.dtype}")
print(f"  min    : {X.min():.3f}")
print(f"  max    : {X.max():.3f}")
print(f"  mean   : {X.mean():.3f}")
if hasattr(X, "toarray"):
    sample = X[:5, :5].toarray()
else:
    sample = X[:5, :5]
print(f"  first 5x5 values:\n{sample}\n")

print("=== .raw ===")
if adata.raw is not None:
    print(f"  shape  : {adata.raw.shape}")
    print(f"  X min  : {adata.raw.X.min():.3f}")
    print(f"  X max  : {adata.raw.X.max():.3f}")
    print(f"  X mean : {adata.raw.X.mean():.3f}")
else:
    print("  NO .raw stored")
print()

print("=== .layers keys ===")
for k in adata.layers.keys():
    L = adata.layers[k]
    print(f"  {k:20s}  min={L.min():.2f}  max={L.max():.2f}  mean={L.mean():.2f}")
if not adata.layers.keys():
    print("  (no layers)")
print()

print("=== obsm keys ===")
for k in adata.obsm.keys():
    print(f"  {k:20s}  shape={adata.obsm[k].shape}")
