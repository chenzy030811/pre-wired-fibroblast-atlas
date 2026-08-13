"""
Re-cluster KC scANVI latent space at lower resolution
to reduce from 18 clusters down to ~13 clusters.

Only re-runs Leiden — does NOT retrain scANVI.
Overwrites leiden_scanvi in place.
"""
import scanpy as sc
from pathlib import Path

path = Path("/Users/apple/Downloads/研究生毕业论文/Cross_Species/"
            "scvi_out/Keratinocyte_scanvi/integrated.h5ad")

print(">>> Loading KC scANVI integrated.h5ad ...")
adata = sc.read_h5ad(path)
adata.obs_names_make_unique()

# neighbors on X_scANVI (should exist already, but recompute to be safe)
print(">>> Rebuilding neighbours on X_scANVI ...")
sc.pp.neighbors(adata, use_rep="X_scANVI")

# try a few resolutions and pick the one closest to 13
target = 13
best_res, best_k, best_diff = None, None, 999

for res in [0.20, 0.25, 0.30, 0.35, 0.40]:
    sc.tl.leiden(adata, resolution=res,
                 key_added=f"_tmp_leiden_{res}",
                 flavor="igraph", n_iterations=2, directed=False)
    k = adata.obs[f"_tmp_leiden_{res}"].nunique()
    print(f"    res={res:.2f}  →  {k} clusters")
    diff = abs(k - target)
    if diff < best_diff:
        best_diff, best_res, best_k = diff, res, k

print(f"\n>>> Best resolution: {best_res}  →  {best_k} clusters (target={target})")

# apply chosen resolution as the canonical leiden_scanvi
adata.obs["leiden_scanvi"] = adata.obs[f"_tmp_leiden_{best_res}"].values

# drop tmp columns
for res in [0.20, 0.25, 0.30, 0.35, 0.40]:
    del adata.obs[f"_tmp_leiden_{res}"]

# save back
adata.write_h5ad(path)
print(f">>> Wrote {path}")
print(f"    Final leiden_scanvi: {adata.obs['leiden_scanvi'].nunique()} clusters")
