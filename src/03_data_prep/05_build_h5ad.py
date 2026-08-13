"""
03b -- Build merged h5ad per compartment from 10x-format exports.

For each compartment (Fibroblast, Keratinocyte):
  - Load 4 species' 10x + metadata
  - Intersect gene set (should already match after ortholog step)
  - Concatenate into one AnnData with 'species' batch column
  - Save as h5ad ready for scVI

Output:
  Cross_Species/h5ad/merged_Fibroblast.h5ad
  Cross_Species/h5ad/merged_Keratinocyte.h5ad
"""
import os
from pathlib import Path
import scipy.io as sio
import scipy.sparse as sp
import pandas as pd
import gzip
import anndata as ad

base = Path("/Users/apple/Downloads/研究生毕业论文/Cross_Species")
tenx_root = base / "tenx_export"
h5ad_dir  = base / "h5ad"
h5ad_dir.mkdir(parents=True, exist_ok=True)

def read_10x(dir_path):
    """Read a 10x-format directory back into an AnnData."""
    dir_path = Path(dir_path)
    with gzip.open(dir_path / "matrix.mtx.gz", "rb") as f:
        mat = sio.mmread(f).tocsr()   # genes x cells

    with gzip.open(dir_path / "barcodes.tsv.gz", "rt") as f:
        barcodes = [l.strip() for l in f]

    with gzip.open(dir_path / "features.tsv.gz", "rt") as f:
        features = [l.split("\t")[0] for l in f]

    meta = pd.read_csv(dir_path / "metadata.csv", index_col=0)
    meta = meta.loc[barcodes]

    adata = ad.AnnData(
        X   = mat.T.tocsr(),  # cells x genes
        obs = meta,
        var = pd.DataFrame(index=features),
    )
    return adata

for comp in ["Fibroblast", "Keratinocyte"]:
    print(f"\n{'='*40}\nCompartment: {comp}\n{'='*40}")

    species_dirs = sorted((tenx_root / comp).iterdir())
    per_species = []

    for sp_dir in species_dirs:
        print(f"  reading {sp_dir.name} ...")
        a = read_10x(sp_dir)
        print(f"    {a.n_obs} cells x {a.n_vars} genes")
        per_species.append(a)

    # concatenate — inner join on genes (should already match)
    merged = ad.concat(per_species, join="inner", label="species_batch",
                       keys=[d.name for d in species_dirs])
    # 'species' column already in each obs from metadata.csv; keep that
    print(f"\n  merged: {merged.n_obs} cells x {merged.n_vars} genes")
    print(merged.obs["species"].value_counts())

    out_path = h5ad_dir / f"merged_{comp}.h5ad"
    merged.write_h5ad(out_path)
    print(f"  wrote {out_path}")

print("\n=== DONE ===")
