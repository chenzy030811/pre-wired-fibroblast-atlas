"""
05d — Marker analysis on Fibroblast v2 scANVI (5-group data, winner).

Input:  scvi_out/Fibroblast_v2_scanvi/integrated.h5ad
Output:
  results/fib_v2_cluster_markers.csv
  preview_output/PREVIEW_v2_Fib_marker_heatmap.png
  preview_output/PREVIEW_v2_Fib_marker_dotplot.png
"""
import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy import sparse

# ============ palette ============
PALETTE = ["#8FB4BE", "#AFC9CF", "#D5E1E3", "#FFFFFF",
           "#EBBFC2", "#E28187", "#D93F49"]
CMAP_PB = LinearSegmentedColormap.from_list("pink_blue", PALETTE, N=256)

plt.rcParams.update({
    "font.family":    "Helvetica", "font.size": 10,
    "axes.linewidth": 0.6,
    "figure.facecolor":"white", "savefig.facecolor":"white",
})

# ============ paths ============
base    = Path("/Users/apple/Downloads/研究生毕业论文/Cross_Species")
in_path = base / "scvi_out/Fibroblast_v2_scanvi/integrated.h5ad"
results = base / "results"; results.mkdir(exist_ok=True, parents=True)
preview = Path("/Users/apple/Downloads/preview_output")

# ============ load + normalize ============
print(">>> Loading Fib v2 scANVI ...")
adata = sc.read_h5ad(in_path)
adata.obs_names_make_unique()
print(f"    {adata.n_obs} cells × {adata.n_vars} HVG")

adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# force numeric ordering
_uniq = sorted(adata.obs["leiden_scanvi"].astype(str).unique(), key=int)
adata.obs["leiden_scanvi"] = pd.Categorical(
    adata.obs["leiden_scanvi"].astype(str),
    categories=_uniq, ordered=True)
print(f"    {len(_uniq)} clusters: {_uniq}")

# ============ rank_genes_groups ============
print(">>> Wilcoxon rank_genes_groups ...")
sc.tl.rank_genes_groups(adata, "leiden_scanvi", method="wilcoxon",
                         pts=True, groups=_uniq)

# collect top 30 per cluster
groups = adata.uns["rank_genes_groups"]["names"].dtype.names
rows = []
for g in groups:
    names  = adata.uns["rank_genes_groups"]["names"][g][:30]
    lfc    = adata.uns["rank_genes_groups"]["logfoldchanges"][g][:30]
    pvals  = adata.uns["rank_genes_groups"]["pvals_adj"][g][:30]
    for rank, (n, l, p) in enumerate(zip(names, lfc, pvals), 1):
        rows.append({"cluster": g, "rank": rank, "gene": n,
                     "log2fc": round(float(l), 3),
                     "pval_adj": float(p)})
df_all = pd.DataFrame(rows)
df_all.to_csv(results / "fib_v2_cluster_markers.csv", index=False)
print(f"    → {results}/fib_v2_cluster_markers.csv")

# ============ heatmap: top 5 per cluster ============
print(">>> Heatmap (top 5 per cluster) ...")
marker_dict = {c: df_all[df_all["cluster"]==c].head(5)["gene"].tolist() for c in _uniq}
sc.pl.matrixplot(
    adata, var_names=marker_dict, groupby="leiden_scanvi",
    categories_order=_uniq, cmap=CMAP_PB,
    standard_scale="var", dendrogram=False, swap_axes=True,
    figsize=(0.55*len(_uniq) + 2.5, 12), show=False,
)
plt.savefig(preview / "PREVIEW_v2_Fib_marker_heatmap.png", dpi=400, bbox_inches="tight")
plt.close()
print(f"    → {preview}/PREVIEW_v2_Fib_marker_heatmap.png")

# ============ dotplot: canonical fib + hero genes ============
canonical = {
    "Universal":     ["PDGFRA", "COL1A1", "DCN", "LUM"],
    "Papillary":     ["DPP4", "WIF1", "COL18A1"],
    "Reticular":     ["DLK1", "TGFBI", "MFAP5", "PI16"],
    "Progenitor":    ["CD34", "PI16", "MFAP5"],
    "Myofibroblast": ["ACTA2", "TAGLN", "POSTN", "CTHRC1"],
    "Adipogenic":    ["PPARG", "LPL"],
    "Hero (Lab)":    ["WNT2", "RSPO3", "GREM1", "SULF1",
                      "CHRDL1", "IL13RA2", "COL13A1"],
}
seen, ordered = set(), []
for cat, genes in canonical.items():
    for g in genes:
        if g in adata.var_names and g not in seen:
            seen.add(g); ordered.append(g)
missing_hero = [g for g in canonical["Hero (Lab)"] if g not in adata.var_names]
print(f"    {len(ordered)} markers present; hero missing: {missing_hero}")

def custom_dotplot(adata, genes, groupby, cats_ordered, cmap, out_path, title):
    gene_idx = [list(adata.var_names).index(g) for g in genes]
    grps = adata.obs[groupby].astype(str).values
    X = adata.X
    if sparse.issparse(X): X = X.toarray()
    X = X[:, gene_idx]
    n_g, n_c = len(cats_ordered), len(genes)
    mean_exp = np.zeros((n_g, n_c)); pct_exp = np.zeros((n_g, n_c))
    for i, grp in enumerate(cats_ordered):
        m = grps == grp
        sub = X[m]
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
    ax.set_xticklabels(cats_ordered, fontsize=10)
    ax.set_yticks(range(n_c)); ax.set_yticklabels(genes[::-1], fontsize=10)
    ax.set_xlim(-0.6, n_g - 0.4); ax.set_ylim(-0.6, n_c - 0.4)
    for s_ in ["top", "right"]: ax.spines[s_].set_visible(False)

    import matplotlib as mpl
    norm = mpl.colors.Normalize(vmin=0, vmax=1)
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
                      cax=axcb, orientation="vertical")
    cb.set_label("Avg.exp.\nscaled", fontsize=10, rotation=0, labelpad=25, va="center")
    cb.ax.tick_params(labelsize=8); cb.outline.set_linewidth(0.5)

    sizes = [20, 40, 60, 80, 100]
    axsz.set_xlim(0,1); axsz.set_ylim(-0.5, len(sizes)-0.3)
    for k, pct in enumerate(sizes):
        s = (pct / 100.0) * max_dot + 5
        axsz.scatter(0.25, k, s=s, c="white", edgecolors="black", linewidths=0.4)
        axsz.text(0.6, k, str(pct), fontsize=9, va="center")
    axsz.text(0.42, len(sizes)-0.05, "Pct.exp", fontsize=10, ha="center")
    axsz.axis("off")

    fig.savefig(out_path, dpi=400, bbox_inches="tight")
    plt.close(fig)

print(">>> Dotplot (canonical + hero) ...")
custom_dotplot(adata, ordered, "leiden_scanvi", _uniq, CMAP_PB,
               preview / "PREVIEW_v2_Fib_marker_dotplot.png",
               "Canonical fibroblast markers × clusters")
print(f"    → {preview}/PREVIEW_v2_Fib_marker_dotplot.png")

# ============ TOP 5 summary ============
print("\n=== TOP 5 MARKERS PER FIB CLUSTER ===")
for g in _uniq:
    top5 = df_all[df_all["cluster"]==g].head(5)["gene"].tolist()
    print(f"  cluster {g:>2}:  " + ", ".join(top5))

print("\n=== DONE ===")
