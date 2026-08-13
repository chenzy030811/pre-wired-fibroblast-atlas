"""
04c -- scVI training on v2 merged data (5 groups: Reindeer_Antler/Back + Human/Mus/Acomys).

Usage:
  python 04c_train_scvi_v2.py Fibroblast
  python 04c_train_scvi_v2.py Keratinocyte

Output:
  scvi_out/<comp>_v2_scvi/integrated.h5ad
  scvi_out/<comp>_v2_scvi/model/
"""
import sys, os
from pathlib import Path
import scanpy as sc
import scvi
import numpy as np

comp = sys.argv[1]
print(f"=== scVI training v2 :: {comp} ===")

root    = Path("/exports/eddie/scratch/s2780961/CrossSpecies")
in_path = root / f"h5ad/merged_{comp}_v2.h5ad"
out_dir = root / f"scvi_out/{comp}_v2_scvi"
out_dir.mkdir(parents=True, exist_ok=True)

print(">>> Loading merged v2 h5ad ...")
adata = sc.read_h5ad(in_path)
adata.obs_names_make_unique()
print(f"    {adata.n_obs} cells x {adata.n_vars} genes")
print(f"    5 groups: {adata.obs['species'].value_counts().to_dict()}")

# ---------- HVG (seurat_v3 needs raw counts) ----------
print("\n>>> Selecting 3000 HVG (seurat_v3, batch=species) ...")
sc.pp.highly_variable_genes(
    adata, n_top_genes=3000,
    flavor="seurat_v3",
    batch_key="species",
    subset=True,
)
print(f"    kept {adata.n_vars} HVG")

# ---------- setup + train ----------
print("\n>>> Setup scVI ...")
scvi.model.SCVI.setup_anndata(adata, batch_key="species")
model = scvi.model.SCVI(adata, n_latent=30)

print(">>> Training scVI (200 epochs) ...")
model.train(max_epochs=200)

# ---------- latent + UMAP + Leiden ----------
print("\n>>> Extracting latent + UMAP + Leiden...")
adata.obsm["X_scVI"] = model.get_latent_representation()
sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5, key_added="leiden_scvi",
             flavor="igraph", n_iterations=2, directed=False)
print(f"    scVI found {adata.obs['leiden_scvi'].nunique()} clusters")

# ---------- save ----------
print("\n>>> Saving ...")
adata.write_h5ad(out_dir / "integrated.h5ad")
model.save(str(out_dir / "model"), save_anndata=False, overwrite=True)
print(f"    Wrote {out_dir/'integrated.h5ad'} + model")
print("\n=== DONE ===")
