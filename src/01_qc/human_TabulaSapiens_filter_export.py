"""
Filter Tabula Sapiens skin h5ad and export to 10x-compatible format.

Filters applied:
  - assay == "10x 3' v3"           (drop Smart-seq2, different chemistry)
  - tissue != "buccal mucosa"      (drop oral mucosa, only keep skin)
  - disease == "normal"            (already all normal, kept for safety)

Output:
  matrix.mtx.gz       - cell x gene sparse counts (10x format)
  barcodes.tsv.gz     - one cell barcode per line
  features.tsv.gz     - gene_id <tab> gene_symbol per line
  metadata.csv        - per-cell metadata (donor, tissue, cell_type, etc.)
"""
import anndata as ad
import scipy.io as sio
import scipy.sparse as sp
import pandas as pd
import gzip
import shutil
from pathlib import Path

raw_dir = Path("/Users/apple/Downloads/研究生毕业论文/human_data/Tabula_Sapiens_raw")
in_h5  = raw_dir / "TabulaSapiens_skin.h5ad"
out_dir = raw_dir
out_dir.mkdir(exist_ok=True)

print(f"Loading {in_h5} ...")
adata = ad.read_h5ad(in_h5)
print(f"  Full dataset: {adata.shape[0]} cells x {adata.shape[1]} genes\n")

print("Applying filters...")
mask = (
    (adata.obs["assay"] == "10x 3' v3") &
    (adata.obs["tissue"] != "buccal mucosa") &
    (adata.obs["disease"] == "normal")
)
adata = adata[mask].copy()
print(f"  After filter: {adata.shape[0]} cells x {adata.shape[1]} genes")
print(f"  Donors        : {adata.obs['donor_id'].nunique()}")
print(f"  Per-donor     :")
print(adata.obs["donor_id"].value_counts())
print(f"  Tissues       :")
print(adata.obs["tissue"].value_counts())
print()

# .X may be raw counts or already log-normalized; check & pick raw layer
if "counts" in adata.layers:
    X = adata.layers["counts"]
    print("Using adata.layers['counts'] (true raw)")
elif "raw_counts" in adata.layers:
    X = adata.layers["raw_counts"]
    print("Using adata.layers['raw_counts']")
elif adata.raw is not None:
    X = adata.raw.X
    print("Using adata.raw.X")
else:
    X = adata.X
    print("Using adata.X (may not be raw -- check carefully)")

# Make sure it's sparse and oriented as cell x gene
if not sp.issparse(X):
    X = sp.csr_matrix(X)

# 10x MTX format expects gene x cell. Transpose to standard 10x convention.
X_t = X.T.tocsr()
print(f"  Matrix orientation: {X_t.shape[0]} genes x {X_t.shape[1]} cells")

print("\nWriting 10x-format files...")

# matrix.mtx (then gzip)
mtx_path = out_dir / "matrix.mtx"
sio.mmwrite(str(mtx_path), X_t, field="integer")
with open(mtx_path, "rb") as f_in, gzip.open(str(mtx_path) + ".gz", "wb") as f_out:
    shutil.copyfileobj(f_in, f_out)
mtx_path.unlink()
print(f"  matrix.mtx.gz       ({(out_dir / 'matrix.mtx.gz').stat().st_size / 1e6:.1f} MB)")

# barcodes.tsv.gz
with gzip.open(out_dir / "barcodes.tsv.gz", "wt") as f:
    for bc in adata.obs_names:
        f.write(bc + "\n")
print(f"  barcodes.tsv.gz     ({len(adata.obs_names):,} cells)")

# features.tsv.gz -- use gene symbols (var index) + gene_id if available
if "feature_name" in adata.var.columns:
    gene_symbols = adata.var["feature_name"].astype(str).tolist()
else:
    gene_symbols = adata.var_names.astype(str).tolist()
gene_ids = adata.var_names.astype(str).tolist()

with gzip.open(out_dir / "features.tsv.gz", "wt") as f:
    for gid, sym in zip(gene_ids, gene_symbols):
        f.write(f"{gid}\t{sym}\tGene Expression\n")
print(f"  features.tsv.gz     ({len(gene_ids):,} genes)")

# metadata.csv
meta_cols = ["donor_id", "tissue", "anatomical_position", "method",
             "cell_type", "broad_cell_class", "compartment",
             "n_genes_by_counts", "total_counts", "pct_counts_mt",
             "sex", "development_stage", "self_reported_ethnicity"]
keep = [c for c in meta_cols if c in adata.obs.columns]
meta = adata.obs[keep].copy()
meta.index.name = "barcode"
meta.to_csv(out_dir / "metadata.csv")
print(f"  metadata.csv        ({meta.shape[0]} rows x {meta.shape[1]} cols)")

print("\nDone.")
