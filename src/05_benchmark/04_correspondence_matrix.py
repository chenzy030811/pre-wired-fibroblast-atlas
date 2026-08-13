"""
08 -- Cross-species correspondence matrix

For each (species × cluster) combination, compute pseudo-bulk mean expression
across all HVG. Then Pearson-correlate every (species, cluster) with every
other (species, cluster).

STRONG correlation between Antler-cluster9 and Human-cluster9 → they are
truly the same fibroblast subtype across species (biological correspondence).

Output:
  preview_output/PREVIEW_v2_correspondence_matrix.png
  results/correspondence_matrix.csv
"""
import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy import sparse

PALETTE = ["#8FB4BE", "#AFC9CF", "#D5E1E3", "#FFFFFF",
           "#EBBFC2", "#E28187", "#D93F49"]
CMAP_PB = LinearSegmentedColormap.from_list("pink_blue", PALETTE, N=256)

SPECIES_ORDER = ["Reindeer_Antler", "Acomys", "Mus", "Human", "Reindeer_Back"]

base    = Path("/Users/apple/Downloads/研究生毕业论文/Cross_Species")
in_path = base / "scvi_out/Fibroblast_v2_scanvi/integrated.h5ad"
results = base / "results"; results.mkdir(exist_ok=True, parents=True)
preview = Path("/Users/apple/Downloads/preview_output")

print(">>> Loading + log-norm ...")
adata = sc.read_h5ad(in_path); adata.obs_names_make_unique()
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sp_arr = adata.obs["species"].astype(str).values
cl_arr = adata.obs["leiden_scanvi"].astype(str).values
_uniq_c = sorted(np.unique(cl_arr), key=int)

X = adata.X
if sparse.issparse(X): X = X.toarray()

# ==== compute pseudo-bulk for each (species, cluster) with >=5 cells ====
print(">>> Computing pseudo-bulks ...")
records, means = [], []
for sp in SPECIES_ORDER:
    for cl in _uniq_c:
        m = (sp_arr == sp) & (cl_arr == cl)
        if m.sum() < 5: continue
        records.append(f"{sp}·{cl}")
        means.append(X[m].mean(axis=0))
pb = np.array(means)  # shape (n_pb, n_genes)
print(f"    {pb.shape[0]} valid (species,cluster) pseudo-bulks")

# ==== Pearson correlation ====
print(">>> Correlating ...")
corr = np.corrcoef(pb)
corr_df = pd.DataFrame(corr, index=records, columns=records)
corr_df.to_csv(results / "correspondence_matrix.csv")

# ==== plot ====
fig, ax = plt.subplots(figsize=(0.4*len(records)+2, 0.4*len(records)+2))
im = ax.imshow(corr, cmap=CMAP_PB, aspect="equal",
               vmin=0.5, vmax=1.0, interpolation="nearest")
ax.set_xticks(range(len(records)))
ax.set_yticks(range(len(records)))
ax.set_xticklabels(records, rotation=90, fontsize=7)
ax.set_yticklabels(records, fontsize=7)
for s_ in ["top","right","left","bottom"]:
    ax.spines[s_].set_visible(False)
ax.tick_params(length=0)

cbar = fig.colorbar(im, ax=ax, shrink=0.5, aspect=15)
cbar.set_label("Pearson r", fontsize=9)
cbar.ax.tick_params(labelsize=8)

plt.title("Cross-species cluster correspondence  "
          "(pseudo-bulk correlation across 5 species × clusters)",
          fontsize=11, fontweight="bold", pad=10)
plt.tight_layout()
plt.savefig(preview / "PREVIEW_v2_correspondence_matrix.png",
            dpi=400, bbox_inches="tight")
plt.close()
print(f"    → {preview}/PREVIEW_v2_correspondence_matrix.png")
print(f"    → {results}/correspondence_matrix.csv")
print("=== DONE ===")
