#!/usr/bin/env python
"""
04 -- Train scVI on merged cross-species h5ad.
Usage:  python 04_train_scvi.py Fibroblast
        python 04_train_scvi.py Keratinocyte
"""
import os
import sys
import scanpy as sc
import anndata as ad
import numpy as np
import scvi
import torch

# --- args ---
if len(sys.argv) != 2 or sys.argv[1] not in ("Fibroblast", "Keratinocyte"):
    sys.exit("Usage: python 04_train_scvi.py {Fibroblast | Keratinocyte}")
comp = sys.argv[1]

# --- paths ---
base = os.environ.get("CROSS_ROOT",
                      "/exports/eddie/scratch/YOURUSER/CrossSpecies")
in_h5   = f"{base}/h5ad/merged_{comp}.h5ad"
out_dir = f"{base}/scvi_out/{comp}_scvi"
os.makedirs(out_dir, exist_ok=True)

print(f">>> Loading {in_h5} ...")
adata = sc.read_h5ad(in_h5)
print(f"    {adata.shape[0]} cells  x  {adata.shape[1]} genes")
print(f"    Species counts:")
print(adata.obs["species"].value_counts())

# --- HVG filtering (top 3000 across all species) ---
print("\n>>> Selecting highly variable genes (top 3,000)...")
sc.pp.highly_variable_genes(
    adata, n_top_genes=3000, flavor="seurat_v3",
    layer=None,  batch_key="species", subset=True
)
print(f"    Kept {adata.n_vars} HVGs")

# --- scVI setup ---
print("\n>>> Setting up scVI...")
scvi.model.SCVI.setup_anndata(adata, batch_key="species")

model = scvi.model.SCVI(
    adata,
    n_layers=2,
    n_latent=30,
    gene_likelihood="nb",
)

# --- Train ---
device = "cuda" if torch.cuda.is_available() else "cpu"
print(f"\n>>> Training on {device}...")
model.train(
    max_epochs=200,
    early_stopping=True,
    early_stopping_patience=20,
    accelerator=device,
)

# --- Save embeddings + model ---
print("\n>>> Extracting latent + running UMAP...")
adata.obsm["X_scVI"] = model.get_latent_representation()
sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=30)
sc.tl.umap(adata, min_dist=0.3)
sc.tl.leiden(adata, resolution=0.5, key_added="leiden_scvi")

# UMAP figure (species + leiden)
import matplotlib.pyplot as plt
plt.rcParams["figure.dpi"] = 150

sc.pl.umap(adata, color=["species", "leiden_scvi"],
           save=f"_{comp}_scvi.png", show=False)

# --- Save outputs ---
print("\n>>> Saving...")
adata.write(os.path.join(out_dir, "integrated.h5ad"))
model.save(out_dir, overwrite=True, save_anndata=False)
print(f"    Wrote {out_dir}/integrated.h5ad + model files")
print("\n=== DONE ===")
