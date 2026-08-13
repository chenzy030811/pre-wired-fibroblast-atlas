"""
08 -- scIB benchmark: score scVI vs scANVI integration quality.

Uses scib-metrics (lightweight, well-maintained fork of scIB).

Metrics computed:
  Batch removal  ── iLISI, kBET, graph connectivity, batch ASW
  Bio conservation ── cLISI, NMI (leiden vs leiden), ARI, isolated labels

Output:
  results/scib_benchmark.csv   — raw metric table
  preview_output/PREVIEW_scib_benchmark.png  — bar chart comparison
"""
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

try:
    from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection
except ImportError:
    print("⚠️  scib-metrics not installed. Install with:")
    print("    pip install scib-metrics")
    raise SystemExit(1)

# ============ paths ============
base    = Path("/Users/apple/Downloads/研究生毕业论文/Cross_Species")
results = base / "results"; results.mkdir(exist_ok=True, parents=True)
preview = Path("/Users/apple/Downloads/preview_output")

scvi_path   = base / "scvi_out/Fibroblast_scvi/integrated.h5ad"
scanvi_path = base / "scvi_out/Fibroblast_scanvi/integrated.h5ad"

# ============ load ============
print(">>> Loading scVI integrated ...")
a_scvi = ad.read_h5ad(scvi_path)
a_scvi.obs_names_make_unique()

print(">>> Loading scANVI integrated ...")
a_scanvi = ad.read_h5ad(scanvi_path)
a_scanvi.obs_names_make_unique()

# ============ align cells (both should have same 21,624 cells) ============
assert a_scvi.n_obs == a_scanvi.n_obs, \
    f"cell counts differ: scVI={a_scvi.n_obs} scANVI={a_scanvi.n_obs}"
print(f"    {a_scvi.n_obs} cells (matched)")

# ============ build combined adata with both embeddings ============
# use scVI's version as the base (has counts, species, leiden)
adata = a_scvi.copy()
adata.obsm["X_scVI"]   = a_scvi.obsm["X_scVI"]
adata.obsm["X_scANVI"] = a_scanvi.obsm["X_scANVI"]

# labels for bio-conservation:
# - use scVI's leiden as a "consensus" proxy label (since we lack shared truth)
adata.obs["label"] = adata.obs["leiden_scvi"].astype("category")
adata.obs["batch"] = adata.obs["species"].astype("category")

# ============ run benchmark ============
print("\n>>> Running scIB benchmark (this takes a few minutes) ...")
bm = Benchmarker(
    adata,
    batch_key="batch",
    label_key="label",
    embedding_obsm_keys=["X_scVI", "X_scANVI"],
    n_jobs=-1,
    bio_conservation_metrics=BioConservation(
        isolated_labels=True,
        nmi_ari_cluster_labels_leiden=True,
        silhouette_label=True,
        clisi_knn=True,
    ),
    batch_correction_metrics=BatchCorrection(
        graph_connectivity=True,
        kbet_per_label=True,
        ilisi_knn=True,
        pcr_comparison=True,
        # silhouette_batch removed in scib-metrics >= 0.5.x
    ),
)
bm.benchmark()

# ============ save results ============
df = bm.get_results(min_max_scale=False)
df.to_csv(results / "scib_benchmark.csv")
print(f"\n    → {results}/scib_benchmark.csv")
print(df)

# ============ bar chart comparison ============
print("\n>>> Drawing bar chart ...")
# transpose so metrics are rows, methods are columns
plot_df = df.T
# drop summary rows if present (Bio conservation, Batch correction, Total)
metric_rows = [r for r in plot_df.index
               if r not in {"Bio conservation", "Batch correction",
                            "Total", "Metric Type"}]
plot_df = plot_df.loc[metric_rows]

fig, ax = plt.subplots(figsize=(10, 5))
xpos = np.arange(len(plot_df))
w = 0.35
scvi_vals   = plot_df["X_scVI"].astype(float).values
scanvi_vals = plot_df["X_scANVI"].astype(float).values

ax.bar(xpos - w/2, scvi_vals,   w, label="scVI",   color="#8AB2D0", edgecolor="black", linewidth=0.4)
ax.bar(xpos + w/2, scanvi_vals, w, label="scANVI", color="#E28A7A", edgecolor="black", linewidth=0.4)

ax.set_xticks(xpos)
ax.set_xticklabels(plot_df.index, rotation=35, ha="right", fontsize=9)
ax.set_ylabel("Score (higher = better)", fontsize=10)
# panel-mode: no title (compose adds it)
ax.legend(frameon=False, fontsize=10)
for s in ["top", "right"]: ax.spines[s].set_visible(False)
ax.grid(axis="y", linestyle=":", linewidth=0.4, color="#D0D0D0")
ax.set_axisbelow(True)

plt.tight_layout()
plt.savefig(preview / "PREVIEW_scib_benchmark.png",
            dpi=400, bbox_inches="tight")
plt.close()
print(f"    → {preview}/PREVIEW_scib_benchmark.png")

# ============ summary ============
print("\n=== SUMMARY ===")
if "Total" in df.index:
    print(f"  Total score — scVI:   {df.loc['Total', 'X_scVI']:.3f}")
    print(f"  Total score — scANVI: {df.loc['Total', 'X_scANVI']:.3f}")
if "Bio conservation" in df.index:
    print(f"  Bio conservation — scVI:   {df.loc['Bio conservation', 'X_scVI']:.3f}")
    print(f"  Bio conservation — scANVI: {df.loc['Bio conservation', 'X_scANVI']:.3f}")
if "Batch correction" in df.index:
    print(f"  Batch correction — scVI:   {df.loc['Batch correction', 'X_scVI']:.3f}")
    print(f"  Batch correction — scANVI: {df.loc['Batch correction', 'X_scANVI']:.3f}")
print("\n=== DONE ===")
