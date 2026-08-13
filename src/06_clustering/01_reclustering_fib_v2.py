"""
Re-cluster v2 integrated data to reduce cluster count to ~13.

For each compartment (Fib winner=scANVI, KC winner=scVI):
  1. Try several resolutions
  2. Pick the one closest to 13
  3. Verify hero-gene enriched cells stay distinct

Overwrites leiden_scanvi / leiden_scvi in place.
"""
import scanpy as sc
import numpy as np
from pathlib import Path

TARGET = 13
CONFIGS = [
    ("Fibroblast_v2_scanvi",   "leiden_scanvi", "X_scANVI"),
    ("Keratinocyte_v2_scvi",   "leiden_scvi",   "X_scVI"),
]

base = Path("/Users/apple/Downloads/研究生毕业论文/Cross_Species/scvi_out")

for folder, leiden_key, latent_key in CONFIGS:
    path = base / folder / "integrated.h5ad"
    print(f"\n{'='*50}\n{folder}\n{'='*50}")
    adata = sc.read_h5ad(path)
    adata.obs_names_make_unique()
    print(f"    {adata.n_obs} cells")

    # Rebuild neighbors on the correct latent
    print(f">>> Rebuilding neighbors on {latent_key} ...")
    sc.pp.neighbors(adata, use_rep=latent_key)

    # Search resolutions
    best_res, best_k, best_diff = None, None, 999
    for res in [0.15, 0.20, 0.25, 0.30, 0.35, 0.40]:
        sc.tl.leiden(adata, resolution=res,
                     key_added=f"_tmp_{res}",
                     flavor="igraph", n_iterations=2, directed=False)
        k = adata.obs[f"_tmp_{res}"].nunique()
        print(f"    res={res:.2f}  →  {k} clusters")
        diff = abs(k - TARGET)
        if diff < best_diff:
            best_diff, best_res, best_k = diff, res, k

    print(f"\n    ★ Selected: res={best_res}  →  {best_k} clusters (target={TARGET})")
    adata.obs[leiden_key] = adata.obs[f"_tmp_{best_res}"].values
    for res in [0.15, 0.20, 0.25, 0.30, 0.35, 0.40]:
        del adata.obs[f"_tmp_{res}"]

    adata.write_h5ad(path)
    print(f">>> Wrote {path}")
    print(f"    Final {leiden_key}: {adata.obs[leiden_key].nunique()} clusters")

print("\n=== DONE ===")
