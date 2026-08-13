"""
Convert Reynolds 2021 (E-MTAB-8142) dense TSV counts -> 10x sparse MTX.
The 3.6 GB TSV would need ~7.5 GB RAM as a dense int matrix in R.
After sparse conversion it shrinks to ~200-400 MB on disk.
"""
import pandas as pd
import scipy.sparse as sp
import scipy.io as sio
import gzip, shutil
from pathlib import Path

raw_dir = Path("/Users/apple/Downloads/研究生毕业论文/human_data/Reynolds_raw")
counts_path = raw_dir / "arrayexpress_counts.txt"

print(f"Reading {counts_path} (3.6 GB)...")
# Streaming read in chunks to avoid memory blowup
chunk_size = 1000  # 1000 genes (rows) at a time
chunks = []
gene_names = []

reader = pd.read_csv(counts_path, sep="\t", chunksize=chunk_size,
                     index_col=0, low_memory=False)
n_chunks = 0
for chunk in reader:
    # Cast value columns to int32 (saves RAM); index column unaffected
    chunks.append(sp.csr_matrix(chunk.astype("int32").values))
    gene_names.extend(chunk.index.tolist())
    n_chunks += 1
    if n_chunks % 5 == 0:
        print(f"  Processed {len(gene_names):>6} genes...")

print(f"  Done reading: {len(gene_names):,} genes")
print("  Stacking sparse matrix...")
mat = sp.vstack(chunks)  # genes x cells
print(f"  Shape: {mat.shape}  (genes x cells)")
print(f"  Non-zero: {mat.nnz:,} ({100*mat.nnz/(mat.shape[0]*mat.shape[1]):.2f}%)")

# Read column names (cell barcodes) from header
with open(counts_path) as f:
    header = f.readline().strip().split("\t")
barcodes = header[1:]  # first col is "Gene"
assert len(barcodes) == mat.shape[1], \
    f"Barcode count mismatch: {len(barcodes)} vs {mat.shape[1]}"

print(f"  {len(barcodes):,} cells\n")

# ---- Write 10x format ----
out_dir = raw_dir
print(f"Writing 10x files to {out_dir}/...")

# matrix.mtx.gz
mtx_path = out_dir / "matrix.mtx"
sio.mmwrite(str(mtx_path), mat, field="integer")
with open(mtx_path, "rb") as f_in, gzip.open(str(mtx_path) + ".gz", "wb") as f_out:
    shutil.copyfileobj(f_in, f_out)
mtx_path.unlink()
print(f"  matrix.mtx.gz   ({(out_dir / 'matrix.mtx.gz').stat().st_size / 1e6:.1f} MB)")

# barcodes.tsv.gz
with gzip.open(out_dir / "barcodes.tsv.gz", "wt") as f:
    for bc in barcodes:
        f.write(bc + "\n")
print(f"  barcodes.tsv.gz ({len(barcodes):,} cells)")

# features.tsv.gz
with gzip.open(out_dir / "features.tsv.gz", "wt") as f:
    for g in gene_names:
        f.write(f"{g}\t{g}\tGene Expression\n")
print(f"  features.tsv.gz ({len(gene_names):,} genes)")

print("\nDone.")
