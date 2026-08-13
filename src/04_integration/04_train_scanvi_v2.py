"""
04c -- scANVI training on v2 merged data (5 groups).

Strategy: first train scVI internally to get Leiden proxy labels,
then convert to scANVI (semi-supervised) using those labels.

Usage:
  python 04c_train_scanvi_v2.py Fibroblast
  python 04c_train_scanvi_v2.py Keratinocyte

Output:
  scvi_out/<comp>_v2_scanvi/integrated.h5ad
  scvi_out/<comp>_v2_scanvi/model/
"""
import sys
from pathlib import Path
import scanpy as sc
import scvi
import numpy as np

comp = sys.argv[1]
print(f"=== scANVI training v2 :: {comp} ===")

root    = Path("/exports/eddie/scratch/s2780961/CrossSpecies")
in_path = root / f"h5ad/merged_{comp}_v2.h5ad"
out_dir = root / f"scvi_out/{comp}_v2_scanvi"
out_dir.mkdir(parents=True, exist_ok=True)

print(">>> Loading merged v2 h5ad ...")
adata = sc.read_h5ad(in_path)
adata.obs_names_make_unique()
print(f"    {adata.n_obs} cells x {adata.n_vars} genes")
print(f"    5 groups: {adata.obs['species'].value_counts().to_dict()}")

# ---------- HVG ----------
print("\n>>> Selecting 3000 HVG (seurat_v3, batch=species) ...")
sc.pp.highly_variable_genes(
    adata, n_top_genes=3000, flavor="seurat_v3",
    batch_key="species", subset=True,
)
print(f"    kept {adata.n_vars} HVG")

# ---------- Step 1: scVI init for proxy labels ----------
print("\n>>> [1/3] scVI init (for proxy labels) ...")
scvi.model.SCVI.setup_anndata(adata, batch_key="species")
scvi_init = scvi.model.SCVI(adata, n_latent=30)
scvi_init.train(max_epochs=200)

adata.obsm["X_scVI_init"] = scvi_init.get_latent_representation()
sc.pp.neighbors(adata, use_rep="X_scVI_init")
sc.tl.leiden(adata, resolution=0.5, key_added="cell_type_proxy",
             flavor="igraph", n_iterations=2, directed=False)
adata.obs["cell_type_proxy"] = adata.obs["cell_type_proxy"].astype(str)
print(f"    proxy labels: {adata.obs['cell_type_proxy'].nunique()} clusters")

# ---------- Step 2: rebuild scVI with labels_key ----------
print("\n>>> [2/3] Rebuild scVI with labels_key ...")
scvi.model.SCVI.setup_anndata(
    adata, batch_key="species", labels_key="cell_type_proxy",
)
scvi_model = scvi.model.SCVI(adata, n_latent=30)
scvi_model.train(max_epochs=200)

# ---------- Step 3: convert to scANVI ----------
print("\n>>> [3/3] Converting to scANVI + training ...")
scanvi_model = scvi.model.SCANVI.from_scvi_model(
    scvi_model,
    unlabeled_category="_UNLABELED_",
    labels_key="cell_type_proxy",
)
scanvi_model.train(max_epochs=100)

# ---------- latent + UMAP + Leiden ----------
print("\n>>> Extracting latent + UMAP + Leiden ...")
adata.obsm["X_scANVI"] = scanvi_model.get_latent_representation()
sc.pp.neighbors(adata, use_rep="X_scANVI")
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5, key_added="leiden_scanvi",
             flavor="igraph", n_iterations=2, directed=False)
print(f"    scANVI found {adata.obs['leiden_scanvi'].nunique()} clusters")

# ---------- save ----------
print("\n>>> Saving ...")
adata.write_h5ad(out_dir / "integrated.h5ad")
scanvi_model.save(str(out_dir / "model"), save_anndata=False, overwrite=True)
print(f"    Wrote {out_dir/'integrated.h5ad'} + model")
print("\n=== DONE ===")
