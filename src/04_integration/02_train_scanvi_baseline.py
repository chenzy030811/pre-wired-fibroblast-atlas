"""
04b -- scANVI training on the cross-species Fibroblast (or KC) data.

Strategy:
  scANVI is a semi-supervised extension of scVI. We use the leiden_scvi
  labels (from the earlier scVI run) as proxy cell-type labels — this is
  the standard practical approach when cross-species unified labels are
  not available (matches Ferreira-lab roadmap Slide 12).

Usage:
  python 04b_train_scanvi.py Fibroblast
  python 04b_train_scanvi.py Keratinocyte

Output:
  scvi_out/<comp>_scanvi/integrated.h5ad
  scvi_out/<comp>_scanvi/model/
"""
import sys
from pathlib import Path
import scanpy as sc
import scvi
import numpy as np

comp = sys.argv[1]
print(f"=== scANVI training :: {comp} ===")

root = Path("/exports/eddie/scratch/s2780961/CrossSpecies")
in_path  = root / f"scvi_out/{comp}_scvi/integrated.h5ad"
out_dir  = root / f"scvi_out/{comp}_scanvi"
out_dir.mkdir(parents=True, exist_ok=True)

# ---------- load: reuse scVI's HVG-filtered raw counts + leiden labels ----
print(">>> Loading scVI integrated.h5ad ...")
adata = sc.read_h5ad(in_path)
adata.obs_names_make_unique()
print(f"    {adata.n_obs} cells x {adata.n_vars} HVG")
print(f"    species : {adata.obs['species'].value_counts().to_dict()}")
print(f"    n_clust : {adata.obs['leiden_scvi'].nunique()}")

# proxy labels for scANVI: leiden clusters found by scVI
adata.obs["cell_type_proxy"] = adata.obs["leiden_scvi"].astype(str)

# ---------- setup + train ----------
print("\n>>> Setup scVI (initial) ...")
scvi.model.SCVI.setup_anndata(
    adata,
    batch_key="species",
    labels_key="cell_type_proxy",
)
scvi_model = scvi.model.SCVI(adata, n_latent=30)
print(">>> Training scVI (200 epochs) ...")
scvi_model.train(max_epochs=200)

print("\n>>> Converting to scANVI ...")
scanvi_model = scvi.model.SCANVI.from_scvi_model(
    scvi_model,
    unlabeled_category="_UNLABELED_",  # unused here; all cells are labelled
    labels_key="cell_type_proxy",
)
print(">>> Training scANVI (100 epochs) ...")
scanvi_model.train(max_epochs=100)

# ---------- latent + UMAP + Leiden ----------
print("\n>>> Extracting latent + running UMAP...")
adata.obsm["X_scANVI"] = scanvi_model.get_latent_representation()
sc.pp.neighbors(adata, use_rep="X_scANVI")
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5, key_added="leiden_scanvi")

n_new = adata.obs["leiden_scanvi"].nunique()
print(f"    scANVI found {n_new} clusters (scVI had {adata.obs['leiden_scvi'].nunique()})")

# ---------- save ----------
print("\n>>> Saving ...")
adata.write_h5ad(out_dir / "integrated.h5ad")
scanvi_model.save(str(out_dir / "model"), save_anndata=False, overwrite=True)
print(f"    Wrote {out_dir/'integrated.h5ad'} + model")
print("\n=== DONE ===")
