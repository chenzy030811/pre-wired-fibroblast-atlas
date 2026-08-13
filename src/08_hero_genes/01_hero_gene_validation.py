"""
06 — Cross-species validation of the 7 hero regenerative genes
     (WNT2 / RSPO3 / GREM1 / SULF1 / CHRDL1 / IL13RA2 / COL13A1)

Tests the Ferreira-lab pre-wired hypothesis:
  → Are hero genes constitutively enriched in regenerating species
    (Reindeer + Acomys) vs scarring species (Mus + Human)?

Outputs 7 figures + 1 CSV:
  1. hero_dotplot_species.png     — dotplot of 7 hero × 12 cluster × 4 species
  2. hero_stacked_violin.png      — stacked violin: 7 hero across 4 species
  3. hero_featureplot_umap.png    — 7 UMAPs, one per hero gene
  4. hero_module_score_umap.png   — 1 UMAP + boxplot of module score
  5. cluster_species_crosstab.png — 12 cluster × 4 species proportion heatmap
  6. hero_boxplot_by_species.png  — box plot of each hero gene split by species
  7. hero_module_stats.txt        — Wilcoxon: regen vs scar module scores

  cluster_species_composition.csv — raw numbers
"""
import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patheffects as mpe
from matplotlib.colors import LinearSegmentedColormap
from scipy import sparse, stats

# ================ paths ================
base    = Path("/Users/apple/Downloads/研究生毕业论文/Cross_Species")
in_path = base / "scvi_out/Fibroblast_scvi/integrated.h5ad"
results = base / "results"; results.mkdir(exist_ok=True, parents=True)
preview = Path("/Users/apple/Downloads/preview_output")

# ================ config ================
HERO_GENES = ["WNT2", "RSPO3", "GREM1", "SULF1",
              "CHRDL1", "IL13RA2", "COL13A1"]

SPECIES_ORDER = ["Reindeer", "Acomys", "Mus", "Human"]     # regen | scar
REGEN_SPECIES = ["Reindeer", "Acomys"]
SCAR_SPECIES  = ["Mus", "Human"]

SPECIES_COLORS = {
    "Reindeer": "#E28A7A",   # unified pastel (matches Fib/KC UMAPs)
    "Acomys":   "#D189A5",
    "Mus":      "#B5C275",
    "Human":    "#8AB2D0",
}

# pink-blue palette from user preference
PALETTE = ["#8FB4BE", "#AFC9CF", "#D5E1E3",
           "#FFFFFF",
           "#EBBFC2", "#E28187", "#D93F49"]
CMAP_PB = LinearSegmentedColormap.from_list("pink_blue", PALETTE, N=256)

plt.rcParams.update({
    "font.family":    "Helvetica",
    "font.size":      10,
    "axes.linewidth": 0.6,
    "figure.facecolor":  "white",
    "savefig.facecolor": "white",
})

# ================ load ================
print(">>> Loading integrated.h5ad ...")
adata = sc.read_h5ad(in_path)
adata.obs_names_make_unique()
print(f"    {adata.n_obs} cells × {adata.n_vars} HVG")

# normalise + log1p (X is raw counts)
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# only hero genes that survived HVG filter
hero_present = [g for g in HERO_GENES if g in adata.var_names]
hero_missing = [g for g in HERO_GENES if g not in adata.var_names]
print(f"    Hero in HVG: {hero_present}")
if hero_missing:
    print(f"    ⚠️  Missing from HVG: {hero_missing}")

# ================ 1. cluster × species crosstab ================
print("\n>>> [1/7] cluster × species crosstab ...")
ct = pd.crosstab(adata.obs["leiden_scvi"].astype(str),
                 adata.obs["species"].astype(str))
ct = ct[[s for s in SPECIES_ORDER if s in ct.columns]]
ct.index = ct.index.astype(int).sort_values().astype(str)
ct = ct.reindex(sorted(ct.index, key=int))
ct.to_csv(results / "cluster_species_composition.csv")
print(f"    → {results}/cluster_species_composition.csv")

# clean heatmap of proportions — matches fib matrixplot style
prop = ct.div(ct.sum(axis=1), axis=0)
# transpose so species are rows (short axis), clusters are columns (long axis)
prop_T = prop.T

fig_w = 0.55 * prop_T.shape[1] + 2.0
fig_h = 0.55 * prop_T.shape[0] + 2.0
fig = plt.figure(figsize=(fig_w, fig_h))
ax   = fig.add_axes([0.14, 0.25, 0.75, 0.58])
axcb = fig.add_axes([0.92, 0.30, 0.015, 0.35])

im = ax.imshow(prop_T.values, cmap=CMAP_PB, aspect="equal",
               vmin=0, vmax=1, interpolation="nearest")
# no gaps between cells

# tick labels
ax.set_xticks(range(prop_T.shape[1]))
ax.set_xticklabels(prop_T.columns, fontsize=10)
ax.set_yticks(range(prop_T.shape[0]))
ax.set_yticklabels(prop_T.index, fontsize=10)
ax.set_xlabel("Fibroblast cluster", fontsize=10, labelpad=6)
ax.set_ylabel("Species", fontsize=10, labelpad=6)
# panel-mode: no internal title (compose adds it)
for s_ in ["top", "right", "left", "bottom"]:
    ax.spines[s_].set_visible(False)
ax.tick_params(length=0)

# annotate % in each cell
for i in range(prop_T.shape[0]):
    for j in range(prop_T.shape[1]):
        v = prop_T.iloc[i, j]
        color = "white" if v > 0.55 else "#222222"
        ax.text(j, i, f"{v*100:.0f}", ha="center", va="center",
                fontsize=8.5, color=color, fontweight="medium")

# vertical colourbar
cbar = fig.colorbar(im, cax=axcb)
cbar.set_label("Proportion", fontsize=9, labelpad=8)
cbar.ax.tick_params(labelsize=8, length=2, width=0.4)
cbar.outline.set_linewidth(0.4)

fig.savefig(preview / "PREVIEW_cluster_species_crosstab.png",
            dpi=400, bbox_inches="tight")
plt.close(fig)
print(f"    → {preview}/PREVIEW_cluster_species_crosstab.png")

# ================ 2. Dotplot: hero × cluster (split by species not possible in one panel;
#                  instead 4 stacked panels, one per species) ================
print("\n>>> [2/7] hero dotplot (per-species panels) ...")

def compute_dotplot_matrices(adata_sub, genes, groupby):
    """Return (mean_scaled, pct) matrices of shape (n_group, n_gene)."""
    genes_here = [g for g in genes if g in adata_sub.var_names]
    idx = [list(adata_sub.var_names).index(g) for g in genes_here]
    grps = adata_sub.obs[groupby].astype(str).values
    unique = sorted(np.unique(grps), key=lambda x: int(x))
    X = adata_sub.X
    if sparse.issparse(X): X = X.toarray()
    X = X[:, idx]
    mean = np.zeros((len(unique), len(genes_here)))
    pct  = np.zeros((len(unique), len(genes_here)))
    for i, g in enumerate(unique):
        m = grps == g
        if m.sum() == 0: continue
        sub = X[m]
        mean[i] = sub.mean(axis=0)
        pct[i]  = (sub > 0).mean(axis=0) * 100
    # scale per gene 0..1 within this species (so colour = relative)
    lo = mean.min(axis=0); hi = mean.max(axis=0)
    rng = np.where(hi - lo == 0, 1, hi - lo)
    scaled = (mean - lo) / rng
    return scaled, pct, unique, genes_here

n_hero = len(hero_present)
fig, axes = plt.subplots(1, 4, figsize=(1.9*n_hero + 2, 5.5), sharey=True)

for ax, sp in zip(axes, SPECIES_ORDER):
    mask = adata.obs["species"].astype(str) == sp
    scaled, pct, groups, genes_h = compute_dotplot_matrices(
        adata[mask], hero_present, "leiden_scvi")

    n_g = len(groups); n_c = len(genes_h)
    for j in range(n_c):
        ax.axhline(j, color="#EEEEEE", lw=0.5, ls=":", zorder=0)
    for i in range(n_g):
        ax.axvline(i, color="#F5F5F5", lw=0.4, ls="-", zorder=0)

    max_dot = 110
    for i in range(n_g):
        for j in range(n_c):
            s = (pct[i, j] / 100.0) * max_dot + 3
            c = CMAP_PB(scaled[i, j])
            ax.scatter(i, n_c - 1 - j, s=s, c=[c],
                       edgecolors="black", linewidths=0.2, zorder=3)

    ax.set_xticks(range(n_g))
    ax.set_xticklabels(groups, fontsize=8)
    ax.set_yticks(range(n_c))
    ax.set_yticklabels(genes_h[::-1], fontsize=9)
    ax.set_xlim(-0.6, n_g - 0.4); ax.set_ylim(-0.6, n_c - 0.4)
    for s_ in ["top", "right"]: ax.spines[s_].set_visible(False)
    ax.set_title(sp, fontsize=11, fontweight="bold",
                 color=SPECIES_COLORS[sp])
    ax.tick_params(length=2, width=0.4)

# panel-mode: no suptitle (compose adds it)
plt.tight_layout()
plt.savefig(preview / "PREVIEW_hero_dotplot_species.png",
            dpi=400, bbox_inches="tight")
plt.close(fig)
print(f"    → {preview}/PREVIEW_hero_dotplot_species.png")

# ================ 3. Stacked violin: hero × species ================
print("\n>>> [3/7] stacked violin ...")
adata_hero = adata[:, hero_present].copy()
adata_hero.obs["species"] = pd.Categorical(
    adata_hero.obs["species"].astype(str),
    categories=SPECIES_ORDER, ordered=True)

sc.pl.stacked_violin(
    adata_hero, hero_present, groupby="species",
    swap_axes=True, show=False,
    figsize=(4.2, 0.65*n_hero + 1.2),
    row_palette=None,
    cmap=CMAP_PB, dendrogram=False,
    standard_scale="var",
)
# panel-mode: no suptitle
plt.savefig(preview / "PREVIEW_hero_stacked_violin.png",
            dpi=400, bbox_inches="tight")
plt.close()
print(f"    → {preview}/PREVIEW_hero_stacked_violin.png")

# ================ 4. UMAP feature plots (7 panels, one per hero gene) ================
print("\n>>> [4/7] UMAP feature plots ...")
umap = adata.obsm["X_umap"]
ncol = 4; nrow = int(np.ceil(n_hero / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(ncol*3.2, nrow*3.0),
                         squeeze=False)

for k, g in enumerate(hero_present):
    ax = axes[k // ncol, k % ncol]
    X = adata[:, g].X
    if sparse.issparse(X): X = X.toarray().flatten()
    else: X = np.asarray(X).flatten()

    # order: low first so high overplot
    order = np.argsort(X)
    sc = ax.scatter(umap[order, 0], umap[order, 1], c=X[order],
                    cmap=CMAP_PB, s=1.5, alpha=0.8,
                    edgecolors="none", rasterized=True,
                    vmin=0, vmax=np.percentile(X[X>0], 98) if (X>0).any() else 1)
    ax.set_xticks([]); ax.set_yticks([])
    for s_ in ax.spines.values(): s_.set_visible(False)
    ax.set_title(g, fontsize=11, fontweight="bold")

# hide empty axes
for k in range(n_hero, nrow*ncol):
    axes[k // ncol, k % ncol].axis("off")

# panel-mode: no suptitle
plt.tight_layout()
plt.savefig(preview / "PREVIEW_hero_featureplot_umap.png",
            dpi=400, bbox_inches="tight")
plt.close(fig)
print(f"    → {preview}/PREVIEW_hero_featureplot_umap.png")

# ================ 5. Module score (7 genes → 1 score per cell) ================
print("\n>>> [5/7] module score ...")
# manual, robust module score: per-cell mean of hero-gene log-norm values
X_hero_full = adata[:, hero_present].X
if sparse.issparse(X_hero_full):
    X_hero_full = X_hero_full.toarray()
adata.obs["Hero_score"] = X_hero_full.mean(axis=1)

# UMAP colored by score
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6),
                                gridspec_kw={"width_ratios": [1.2, 1]})

# left: UMAP
score = adata.obs["Hero_score"].values
order = np.argsort(score)
vmin, vmax = np.percentile(score, [1, 99])
sc_ = ax1.scatter(umap[order, 0], umap[order, 1], c=score[order],
                  cmap=CMAP_PB, s=1.5, alpha=0.8,
                  edgecolors="none", rasterized=True,
                  vmin=vmin, vmax=vmax)
ax1.set_xticks([]); ax1.set_yticks([])
for s_ in ax1.spines.values(): s_.set_visible(False)
ax1.set_title("Hero-gene module score (UMAP)", fontsize=11, fontweight="bold")
cb = fig.colorbar(sc_, ax=ax1, shrink=0.6, label="Module score")
cb.outline.set_linewidth(0.4)

# right: boxplot by species
data = [adata.obs.loc[adata.obs["species"]==sp, "Hero_score"].values
        for sp in SPECIES_ORDER]
bp = ax2.boxplot(data, labels=SPECIES_ORDER, patch_artist=True,
                 widths=0.6, showfliers=False,
                 medianprops=dict(color="black", linewidth=1.2))
for patch, sp in zip(bp["boxes"], SPECIES_ORDER):
    patch.set_facecolor(SPECIES_COLORS[sp])
    patch.set_edgecolor("black"); patch.set_linewidth(0.5); patch.set_alpha(0.85)
for element in ["whiskers", "caps"]:
    for line in bp[element]:
        line.set_color("black"); line.set_linewidth(0.6)
ax2.set_ylabel("Hero-gene module score")
ax2.set_title("Module score across species", fontsize=11, fontweight="bold")
ax2.spines["top"].set_visible(False); ax2.spines["right"].set_visible(False)
ax2.grid(axis="y", linestyle=":", linewidth=0.4, color="#D0D0D0")
ax2.set_axisbelow(True)

# --- ADD: statistical annotation (Regen vs Scar Mann-Whitney U) ---
score_arr = adata.obs["Hero_score"].values
sp_arr_all = adata.obs["species"].astype(str).values
regen_scores = score_arr[np.isin(sp_arr_all, REGEN_SPECIES)]
scar_scores  = score_arr[np.isin(sp_arr_all, SCAR_SPECIES)]
u_stat, p_val = stats.mannwhitneyu(regen_scores, scar_scores, alternative="two-sided")

# significance stars
if p_val < 1e-4:   stars = "****"
elif p_val < 1e-3: stars = "***"
elif p_val < 1e-2: stars = "**"
elif p_val < 5e-2: stars = "*"
else:              stars = "ns"

# draw bracket over regen (positions 1,2) vs scar (positions 3,4)
ymax_bp = max(np.percentile(d, 90) for d in data)
y0 = ymax_bp * 1.15
h  = ymax_bp * 0.04
ax2.plot([1, 1, 4, 4], [y0, y0+h, y0+h, y0],
         color="black", lw=1.0, clip_on=False)
ax2.text(2.5, y0 + h*1.5,
         f"{stars}   p = {p_val:.1e}",
         ha="center", va="bottom", fontsize=9, fontweight="bold")

# labels under species: which is regen / scar
ax2.text(1.5, -ymax_bp*0.15, "regenerating",
         ha="center", fontsize=8, color="#B85049", style="italic")
ax2.text(3.5, -ymax_bp*0.15, "scarring",
         ha="center", fontsize=8, color="#4A7CB0", style="italic")

# extend ylim to make room for bracket
ax2.set_ylim(-ymax_bp*0.22, y0 + h*3)

plt.tight_layout()
plt.savefig(preview / "PREVIEW_hero_module_score_umap.png",
            dpi=400, bbox_inches="tight")
plt.close(fig)
print(f"    → {preview}/PREVIEW_hero_module_score_umap.png")

# ================ 6. Box plot: each hero gene split by 4 species ================
print("\n>>> [6/7] per-hero-gene boxplot by species ...")
X_hero = adata[:, hero_present].X
if sparse.issparse(X_hero): X_hero = X_hero.toarray()
sp_arr = adata.obs["species"].astype(str).values

fig, axes = plt.subplots(1, n_hero, figsize=(1.8*n_hero, 3.5), sharey=False)
if n_hero == 1: axes = [axes]

for k, g in enumerate(hero_present):
    ax = axes[k]
    expr = X_hero[:, k]
    data = [expr[sp_arr == sp] for sp in SPECIES_ORDER]
    bp = ax.boxplot(data, labels=SPECIES_ORDER, patch_artist=True,
                    widths=0.6, showfliers=False,
                    medianprops=dict(color="black", linewidth=1.0))
    for patch, sp in zip(bp["boxes"], SPECIES_ORDER):
        patch.set_facecolor(SPECIES_COLORS[sp])
        patch.set_edgecolor("black"); patch.set_linewidth(0.4); patch.set_alpha(0.85)
    for element in ["whiskers", "caps"]:
        for line in bp[element]:
            line.set_color("black"); line.set_linewidth(0.5)
    ax.set_title(g, fontsize=10, fontweight="bold")
    ax.set_xticklabels(SPECIES_ORDER, rotation=35, ha="right", fontsize=8)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=8)
    if k == 0: ax.set_ylabel("log-normalised expression", fontsize=9)

    # --- FIX: force visible y-axis range for sparse genes + annotate % positive ---
    gene_max = expr.max()
    if gene_max < 0.5:
        # sparsely expressed — force ylim so panel isn't visually empty
        ax.set_ylim(-0.05, 1.0)
        pct_pos = (expr > 0).mean() * 100
        ax.text(0.5, 0.92,
                f"{pct_pos:.1f}% cells > 0",
                transform=ax.transAxes, ha="center", va="top",
                fontsize=7, color="#555555", style="italic")
    else:
        ax.set_ylim(-0.05, gene_max * 1.10)

# panel-mode: no suptitle
plt.tight_layout()
plt.savefig(preview / "PREVIEW_hero_boxplot_by_species.png",
            dpi=400, bbox_inches="tight")
plt.close(fig)
print(f"    → {preview}/PREVIEW_hero_boxplot_by_species.png")

# ================ 7. Statistical tests ================
print("\n>>> [7/7] statistical tests (regen vs scar) ...")
lines = ["=== HERO GENE STATISTICS ===", ""]

# module score test
score = adata.obs["Hero_score"].values
sp    = adata.obs["species"].astype(str).values
regen = score[np.isin(sp, REGEN_SPECIES)]
scar  = score[np.isin(sp, SCAR_SPECIES)]
u, p  = stats.mannwhitneyu(regen, scar, alternative="greater")
lines += [
    "Module-score comparison (all fibroblasts):",
    f"  Regenerators (Reindeer+Acomys) n={len(regen)}  median={np.median(regen):+.3f}",
    f"  Scarrers     (Mus+Human)       n={len(scar)}   median={np.median(scar):+.3f}",
    f"  Mann-Whitney U = {u:.1f}",
    f"  p-value (one-sided, regen>scar) = {p:.3e}",
    "",
]

# per-gene test
lines.append("Per-gene comparison (Mann-Whitney, regen > scar):")
lines.append(f"{'gene':10s} {'median_regen':>12s} {'median_scar':>12s} "
             f"{'log2FC':>8s} {'p_value':>10s}")
X_hero = adata[:, hero_present].X
if sparse.issparse(X_hero): X_hero = X_hero.toarray()
for k, g in enumerate(hero_present):
    e = X_hero[:, k]
    r = e[np.isin(sp, REGEN_SPECIES)]
    s = e[np.isin(sp, SCAR_SPECIES)]
    u, p = stats.mannwhitneyu(r, s, alternative="greater")
    lfc = np.log2((r.mean() + 1e-3) / (s.mean() + 1e-3))
    lines.append(f"{g:10s} {np.median(r):+12.3f} {np.median(s):+12.3f} "
                 f"{lfc:+8.3f} {p:10.2e}")

# per-cluster composition summary
lines += ["", "=== KEY OBSERVATIONS FROM CROSSTAB ==="]
for c in ct.index:
    row = ct.loc[c]
    total = row.sum()
    if total == 0: continue
    prop_regen = (row.get("Reindeer",0) + row.get("Acomys",0)) / total
    if prop_regen >= 0.6:
        lines.append(
            f"  cluster {c}: {prop_regen*100:.0f}% regenerator "
            f"(R={row.get('Reindeer',0)}, A={row.get('Acomys',0)}, "
            f"M={row.get('Mus',0)}, H={row.get('Human',0)})"
        )

report = "\n".join(lines)
(results / "hero_module_stats.txt").write_text(report)
print("\n" + report)
print(f"\n    → {results}/hero_module_stats.txt")

print("\n=== DONE ===")
