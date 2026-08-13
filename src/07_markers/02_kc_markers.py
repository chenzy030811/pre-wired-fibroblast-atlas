"""
05c — Marker gene discovery for cross-species Keratinocyte clusters
      (parallel to 05b for Fibroblast).

Input:
  scvi_out/Keratinocyte_scvi/integrated.h5ad  (raw counts, HVG, leiden_scvi)

Steps:
  1. normalize + log1p
  2. rank_genes_groups (Wilcoxon) per cluster
  3. Export CSV of top 30 markers per cluster
  4. Heatmap (top 5 markers per cluster)
  5. Dotplot: canonical keratinocyte subtype markers × N clusters

Output:
  results/kc_cluster_markers.csv
  preview_output/PREVIEW_KC_marker_heatmap.png
  preview_output/PREVIEW_KC_marker_dotplot.png
"""
import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy import sparse

# ============ palette — same pink-to-blue-grey as Fib ==
PALETTE_HEX = ["#8FB4BE", "#AFC9CF", "#D5E1E3",
               "#FFFFFF",
               "#EBBFC2", "#E28187", "#D93F49"]
CMAP_PINK_BLUE = LinearSegmentedColormap.from_list("pink_blue", PALETTE_HEX, N=256)

plt.rcParams.update({
    "font.family":    "Helvetica",
    "font.size":      10,
    "axes.linewidth": 0.6,
    "figure.facecolor":  "white",
    "savefig.facecolor": "white",
})

# ============ paths ============
base    = Path("/Users/apple/Downloads/研究生毕业论文/Cross_Species")
in_path = base / "scvi_out/Keratinocyte_scvi/integrated.h5ad"
results = base / "results"; results.mkdir(exist_ok=True, parents=True)
preview = Path("/Users/apple/Downloads/preview_output")

# ============ load + prep ============
print(">>> Loading Keratinocyte integrated.h5ad ...")
adata = sc.read_h5ad(in_path)
adata.obs_names_make_unique()
print(f"    {adata.n_obs} cells × {adata.n_vars} HVG")

print(">>> Normalising + log1p ...")
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# force numeric ordering of leiden_scvi (0,1,2,...12) instead of string
# ("0","1","10","11","12","2","3",...)
_uniq = sorted(adata.obs["leiden_scvi"].astype(str).unique(), key=int)
adata.obs["leiden_scvi"] = pd.Categorical(
    adata.obs["leiden_scvi"].astype(str),
    categories=_uniq, ordered=True,
)
print(f"    cluster order enforced: {_uniq}")

# ============ marker discovery ============
print(">>> Running rank_genes_groups (Wilcoxon) ...")
sc.tl.rank_genes_groups(
    adata, "leiden_scvi",
    method="wilcoxon", pts=True,
    groups=_uniq,          # explicit numeric order 0,1,2,...,12
)

# --------- CSV ---------
print(">>> Building CSV ...")
groups = adata.uns["rank_genes_groups"]["names"].dtype.names
rows = []
for g in groups:
    names   = adata.uns["rank_genes_groups"]["names"][g][:30]
    scores  = adata.uns["rank_genes_groups"]["scores"][g][:30]
    lfc     = adata.uns["rank_genes_groups"]["logfoldchanges"][g][:30]
    pvals   = adata.uns["rank_genes_groups"]["pvals_adj"][g][:30]
    for rank, (n, s, l, p) in enumerate(zip(names, scores, lfc, pvals), 1):
        rows.append({
            "cluster": g, "rank": rank, "gene": n,
            "score": round(float(s), 3),
            "log2fc": round(float(l), 3),
            "pval_adj": float(p),
        })
csv_path = results / "kc_cluster_markers.csv"
pd.DataFrame(rows).to_csv(csv_path, index=False)
print(f"    → {csv_path}")

# --------- Heatmap: use sc.pl.matrixplot with explicit ordering
#           (rank_genes_groups_matrixplot ignores group order — buggy) ---
print(">>> Drawing marker matrix-heatmap ...")
n_clust = len(_uniq)

# build ordered {cluster -> [top 5 marker genes]}  in numeric 0..12 order
df_all = pd.DataFrame(rows)
marker_dict = {}
for c in _uniq:
    top5 = df_all[df_all["cluster"] == c].head(5)["gene"].tolist()
    marker_dict[c] = top5

sc.pl.matrixplot(
    adata, var_names=marker_dict, groupby="leiden_scvi",
    categories_order=_uniq,         # ← this one actually works
    cmap=CMAP_PINK_BLUE,
    standard_scale="var",
    dendrogram=False, swap_axes=True,
    figsize=(0.55*n_clust + 2.5, 12),
    show=False,
)
plt.suptitle("Top 5 markers per keratinocyte cluster",
             y=1.005, fontsize=13, fontweight="bold")
plt.savefig(preview / "PREVIEW_KC_marker_heatmap.png",
            dpi=400, bbox_inches="tight")
plt.close()
print(f"    → {preview}/PREVIEW_KC_marker_heatmap.png")

# --------- Canonical KC subtype markers ---------
# Drawn from Reynolds 2021, Cheng 2018, Joost 2016, Sole-Boldo 2020,
# and standard skin scRNA-seq references.
canonical_kc = {
    "Basal":         ["KRT5", "KRT14", "TP63", "ITGA6", "ITGB1"],
    "Suprabasal":    ["KRT1", "KRT10", "IVL"],
    "Granular":      ["LOR", "FLG", "KRT2"],
    "Hair follicle": ["KRT15", "KRT17", "SOX9", "LGR5", "LHX2"],
    "Sebaceous":     ["FASN", "SCD", "MGST1"],
    "Cycling":       ["MKI67", "TOP2A", "PCNA", "CDK1"],
    "Sweat":         ["DCD", "MUCL1"],
}
seen, ordered = set(), []
for cat, genes in canonical_kc.items():
    for g in genes:
        if g in adata.var_names and g not in seen:
            seen.add(g); ordered.append(g)

missing = [g for cat, gs in canonical_kc.items() for g in gs
           if g not in adata.var_names]
print(f"    {len(ordered)} canonical KC markers found in HVG")
if missing:
    print(f"    ⚠️  Not in HVG: {missing[:10]}...")

# --------- Custom reference-style dotplot ---------
def custom_dotplot(adata, genes, groupby, cmap, out_path, title):
    genes_present = [g for g in genes if g in adata.var_names]
    gene_idx = [list(adata.var_names).index(g) for g in genes_present]
    groups = adata.obs[groupby].astype(str).values
    unique_groups = sorted(np.unique(groups), key=lambda x: int(x))

    Xall = adata.X
    if sparse.issparse(Xall): Xall = Xall.toarray()
    Xall = Xall[:, gene_idx]

    n_g, n_c = len(unique_groups), len(genes_present)
    mean_exp = np.zeros((n_g, n_c))
    pct_exp  = np.zeros((n_g, n_c))
    for i, grp in enumerate(unique_groups):
        m = groups == grp
        sub = Xall[m]
        mean_exp[i] = sub.mean(axis=0)
        pct_exp[i]  = (sub > 0).mean(axis=0) * 100

    lo = mean_exp.min(axis=0); hi = mean_exp.max(axis=0)
    rng = np.where(hi - lo == 0, 1, hi - lo)
    scaled = (mean_exp - lo) / rng

    fig_w = 0.32*n_g + 3.0
    fig_h = 0.22*n_c + 1.3
    fig = plt.figure(figsize=(fig_w, fig_h))
    ax   = fig.add_axes([0.13, 0.08, 0.63, 0.85])
    axcb = fig.add_axes([0.82, 0.55, 0.018, 0.30])
    axsz = fig.add_axes([0.82, 0.15, 0.10, 0.30])

    for j in range(n_c):
        ax.axhline(j, color="#D9D9D9", lw=0.5, ls=":", zorder=0)
    for i in range(n_g):
        ax.axvline(i, color="#F0F0F0", lw=0.4, ls="-", zorder=0)

    max_dot = 110
    for i in range(n_g):
        for j in range(n_c):
            s = (pct_exp[i, j] / 100.0) * max_dot + 5
            c = cmap(scaled[i, j])
            ax.scatter(i, n_c - 1 - j, s=s, c=[c],
                       edgecolors="black", linewidths=0.25, zorder=3)

    ax.set_xticks(range(n_g))
    ax.set_xticklabels(unique_groups, fontsize=10)
    ax.set_yticks(range(n_c))
    ax.set_yticklabels(genes_present[::-1], fontsize=10)
    ax.set_xlim(-0.6, n_g - 0.4)
    ax.set_ylim(-0.6, n_c - 0.4)
    for s in ["top", "right"]: ax.spines[s].set_visible(False)
    ax.spines["left"].set_linewidth(0.6)
    ax.spines["bottom"].set_linewidth(0.6)
    ax.tick_params(length=2, width=0.5)

    import matplotlib as mpl
    norm = mpl.colors.Normalize(vmin=0, vmax=1)
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
                      cax=axcb, orientation="vertical")
    cb.set_label("Avg.exp.\nscaled", fontsize=10, rotation=0, labelpad=25, va="center")
    cb.ax.tick_params(labelsize=8, length=2, width=0.5)
    cb.outline.set_linewidth(0.5)

    sizes = [20, 40, 60, 80, 100]
    axsz.set_xlim(0, 1); axsz.set_ylim(-0.5, len(sizes) - 0.3)
    for k, pct in enumerate(sizes):
        s = (pct / 100.0) * max_dot + 5
        axsz.scatter(0.25, k, s=s, c="white",
                     edgecolors="black", linewidths=0.4)
        axsz.text(0.6, k, str(pct), fontsize=9, va="center")
    axsz.text(0.42, len(sizes) - 0.05, "Pct.exp",
              fontsize=10, ha="center", fontweight="normal")
    axsz.axis("off")

    fig.suptitle(title, y=0.99, fontsize=13, fontweight="bold")
    fig.savefig(out_path, dpi=400, bbox_inches="tight")
    plt.close(fig)

print(">>> Drawing custom dotplot ...")
custom_dotplot(
    adata, ordered, "leiden_scvi", CMAP_PINK_BLUE,
    preview / "PREVIEW_KC_marker_dotplot.png",
    f"Canonical keratinocyte subtype markers × {n_clust} clusters",
)
print(f"    → {preview}/PREVIEW_KC_marker_dotplot.png")

# --------- Summary ---------
print("\n=== TOP 5 MARKERS PER KC CLUSTER ===")
df = pd.DataFrame(rows)
for g in groups:
    top5 = df[df["cluster"] == g].head(5)["gene"].tolist()
    print(f"  cluster {g:>2}:  " + ", ".join(top5))

print("\n=== DONE ===")
