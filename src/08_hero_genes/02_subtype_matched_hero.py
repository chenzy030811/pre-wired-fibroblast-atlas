"""
07 -- Subtype-matched hero-gene validation on Fib v2 scANVI (winner).

THE actual hypothesis test:
  Within each hero-gene-enriched fibroblast subtype,
  do regenerating species (Reindeer_Antler + Acomys) express
  hero genes at higher level than scarring species
  (Reindeer_Back + Mus + Human)?

This solves the apples-vs-oranges problem of bulk comparison.

Outputs (all in preview_output/ and results/):
  1. PREVIEW_v2_cluster_species_crosstab.png   (14 clusters × 5 groups)
  2. PREVIEW_v2_hero_by_cluster_species.png    (hero × cluster × species heatmap)
  3. PREVIEW_v2_module_score_matched.png       (per-cluster boxplot with p-values)
  4. results/subtype_matched_stats.csv
"""
import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
from scipy import sparse, stats

# ============ config ============
HERO_GENES = ["WNT2", "RSPO3", "GREM1", "SULF1",
              "CHRDL1", "IL13RA2", "COL13A1"]

# 5 groups in canonical order
SPECIES_ORDER = ["Reindeer_Antler", "Acomys", "Mus", "Human", "Reindeer_Back"]
REGEN_SPECIES = ["Reindeer_Antler", "Acomys"]
SCAR_SPECIES  = ["Mus", "Human", "Reindeer_Back"]

SPECIES_COLORS = {
    "Reindeer_Antler": "#E28A7A",
    "Acomys":          "#D189A5",
    "Mus":             "#B5C275",
    "Human":           "#8AB2D0",
    "Reindeer_Back":   "#8A5A55",
}

# palette
PALETTE = ["#8FB4BE", "#AFC9CF", "#D5E1E3", "#FFFFFF",
           "#EBBFC2", "#E28187", "#D93F49"]
CMAP_PB = LinearSegmentedColormap.from_list("pink_blue", PALETTE, N=256)

plt.rcParams.update({
    "font.family":"Helvetica", "font.size":10,
    "axes.linewidth":0.6,
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
print(f"    {adata.n_obs} cells")

adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# force numeric cluster order
_uniq_c = sorted(adata.obs["leiden_scanvi"].astype(str).unique(), key=int)
adata.obs["leiden_scanvi"] = pd.Categorical(
    adata.obs["leiden_scanvi"].astype(str),
    categories=_uniq_c, ordered=True)
print(f"    {len(_uniq_c)} clusters, {len(SPECIES_ORDER)} species groups")

# only hero genes actually present
hero_present = [g for g in HERO_GENES if g in adata.var_names]
missing = [g for g in HERO_GENES if g not in adata.var_names]
print(f"    Hero present: {hero_present}")
if missing:
    print(f"    ⚠️  Missing: {missing}")

# ==============================================================
# [1] Cluster × Species crosstab (heatmap of proportions)
# ==============================================================
print("\n>>> [1] Cluster × Species crosstab ...")
ct = pd.crosstab(adata.obs["leiden_scanvi"].astype(str),
                 adata.obs["species"].astype(str))
ct = ct[[s for s in SPECIES_ORDER if s in ct.columns]]
ct = ct.reindex(_uniq_c)
prop = ct.div(ct.sum(axis=1), axis=0)
prop.to_csv(results / "cluster_species_composition_v2.csv")

# heatmap — transposed so species=rows, clusters=cols
prop_T = prop.T
fig_w = 0.55 * prop_T.shape[1] + 2.0
fig_h = 0.55 * prop_T.shape[0] + 2.0
fig = plt.figure(figsize=(fig_w, fig_h))
ax   = fig.add_axes([0.14, 0.25, 0.75, 0.58])
axcb = fig.add_axes([0.92, 0.30, 0.015, 0.35])

im = ax.imshow(prop_T.values, cmap=CMAP_PB, aspect="equal",
               vmin=0, vmax=1, interpolation="nearest")
ax.set_xticks(range(prop_T.shape[1]))
ax.set_xticklabels(prop_T.columns, fontsize=10)
ax.set_yticks(range(prop_T.shape[0]))
ax.set_yticklabels(prop_T.index, fontsize=10)
ax.set_xlabel("Fibroblast cluster", fontsize=10, labelpad=6)
ax.set_ylabel("Species group", fontsize=10, labelpad=6)
for s_ in ["top","right","left","bottom"]: ax.spines[s_].set_visible(False)
ax.tick_params(length=0)

for i in range(prop_T.shape[0]):
    for j in range(prop_T.shape[1]):
        v = prop_T.iloc[i, j]
        color = "white" if v > 0.55 else "#222222"
        ax.text(j, i, f"{v*100:.0f}", ha="center", va="center",
                fontsize=8.5, color=color, fontweight="medium")

cbar = fig.colorbar(im, cax=axcb)
cbar.set_label("Proportion", fontsize=9, labelpad=8)
cbar.ax.tick_params(labelsize=8); cbar.outline.set_linewidth(0.4)

fig.savefig(preview / "PREVIEW_v2_cluster_species_crosstab.png",
            dpi=400, bbox_inches="tight")
plt.close(fig)
print(f"    → {preview}/PREVIEW_v2_cluster_species_crosstab.png")

# ==============================================================
# [2] Module score per cell, then per-cluster comparison
# ==============================================================
print("\n>>> [2] Computing hero-gene module score per cell ...")
X_hero = adata[:, hero_present].X
if sparse.issparse(X_hero):
    X_hero = X_hero.toarray()
adata.obs["Hero_score"] = X_hero.mean(axis=1)

# ==============================================================
# [3] Subtype-matched stats: for each cluster × each hero-gene,
#     Mann-Whitney U (regen vs scar). Also for module score.
# ==============================================================
print("\n>>> [3] Subtype-matched statistics ...")
sp_arr = adata.obs["species"].astype(str).values
cl_arr = adata.obs["leiden_scanvi"].astype(str).values

rows = []
for cl in _uniq_c:
    m_cl = cl_arr == cl
    if m_cl.sum() < 20: continue
    m_regen = m_cl & np.isin(sp_arr, REGEN_SPECIES)
    m_scar  = m_cl & np.isin(sp_arr, SCAR_SPECIES)
    if m_regen.sum() < 5 or m_scar.sum() < 5:
        continue

    # Hero_score
    r = adata.obs.loc[m_regen, "Hero_score"].values
    s = adata.obs.loc[m_scar,  "Hero_score"].values
    if r.max() == 0 and s.max() == 0:
        u, p = np.nan, np.nan
    else:
        u, p = stats.mannwhitneyu(r, s, alternative="two-sided")
    rows.append({
        "cluster": cl, "gene": "Hero_module",
        "n_regen": int(m_regen.sum()), "n_scar": int(m_scar.sum()),
        "median_regen": float(np.median(r)),
        "median_scar":  float(np.median(s)),
        "log2fc_mean":  float(np.log2((r.mean()+1e-3)/(s.mean()+1e-3))),
        "p_value":      float(p),
    })

    # per hero gene
    for k, g in enumerate(hero_present):
        e_r = X_hero[m_regen, k]
        e_s = X_hero[m_scar, k]
        if e_r.max() == 0 and e_s.max() == 0:
            u, p = np.nan, np.nan
        else:
            u, p = stats.mannwhitneyu(e_r, e_s, alternative="two-sided")
        rows.append({
            "cluster": cl, "gene": g,
            "n_regen": int(m_regen.sum()), "n_scar": int(m_scar.sum()),
            "median_regen": float(np.median(e_r)),
            "median_scar":  float(np.median(e_s)),
            "log2fc_mean":  float(np.log2((e_r.mean()+1e-3)/(e_s.mean()+1e-3))),
            "p_value":      float(p),
        })

stats_df = pd.DataFrame(rows)
stats_df.to_csv(results / "subtype_matched_stats.csv", index=False)
print(f"    → {results}/subtype_matched_stats.csv")

# ==============================================================
# [4] Master boxplot: module score per cluster, split by 5 species
# ==============================================================
print("\n>>> [4] Master boxplot (module score per cluster × 5 species) ...")
# focus on clusters with >= 20 cells overall
valid_clusters = [c for c in _uniq_c if (cl_arr == c).sum() >= 20]

ncol = 5
nrow = int(np.ceil(len(valid_clusters) / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(ncol*2.8, nrow*2.6),
                         sharey=True, squeeze=False)

score_arr = adata.obs["Hero_score"].values

for k, cl in enumerate(valid_clusters):
    r, c = k // ncol, k % ncol
    ax = axes[r, c]
    m_cl = cl_arr == cl
    data = []
    ticks = []
    colors = []
    for sp in SPECIES_ORDER:
        m = m_cl & (sp_arr == sp)
        if m.sum() == 0:
            continue
        data.append(score_arr[m])
        ticks.append(sp)
        colors.append(SPECIES_COLORS[sp])
    if not data:
        ax.axis("off"); continue

    bp = ax.boxplot(data, patch_artist=True, widths=0.6, showfliers=False,
                    medianprops=dict(color="black", linewidth=1.0))
    for patch, col in zip(bp["boxes"], colors):
        patch.set_facecolor(col); patch.set_edgecolor("black")
        patch.set_linewidth(0.4); patch.set_alpha(0.85)
    for element in ["whiskers", "caps"]:
        for line in bp[element]:
            line.set_color("black"); line.set_linewidth(0.5)

    ax.set_xticks(range(1, len(ticks)+1))
    ax.set_xticklabels(ticks, rotation=45, ha="right", fontsize=7)

    # add p-value bracket if module row exists
    row = stats_df[(stats_df["cluster"]==cl) & (stats_df["gene"]=="Hero_module")]
    if not row.empty and not np.isnan(row["p_value"].iloc[0]):
        p_val = row["p_value"].iloc[0]
        lfc   = row["log2fc_mean"].iloc[0]
        if p_val < 1e-4: stars = "****"
        elif p_val < 1e-3: stars = "***"
        elif p_val < 1e-2: stars = "**"
        elif p_val < 5e-2: stars = "*"
        else:              stars = "ns"
        # color by direction (regen>scar is "supports hypothesis")
        c_txt = "#B85049" if lfc > 0 else "#7A7A7A"
        ax.text(0.98, 0.95, f"{stars}\np={p_val:.1e}\nlog2FC={lfc:+.2f}",
                transform=ax.transAxes, ha="right", va="top",
                fontsize=7, color=c_txt, fontweight="bold")

    # cell counts
    n = int(m_cl.sum())
    ax.set_title(f"cluster {cl} · n={n}",
                 fontsize=10, fontweight="bold", pad=4)
    for s_ in ["top","right"]: ax.spines[s_].set_visible(False)
    if c == 0: ax.set_ylabel("Hero-gene\nmodule score", fontsize=9)
    ax.tick_params(labelsize=7)

# hide empty axes
for k in range(len(valid_clusters), nrow*ncol):
    axes[k//ncol, k%ncol].axis("off")

fig.suptitle("Subtype-matched hero-gene module score · regen vs scar test",
             y=1.005, fontsize=13, fontweight="bold")
plt.tight_layout()
plt.savefig(preview / "PREVIEW_v2_module_score_matched.png",
            dpi=400, bbox_inches="tight")
plt.close()
print(f"    → {preview}/PREVIEW_v2_module_score_matched.png")

# ==============================================================
# [5] Print top significant findings
# ==============================================================
print("\n=== TOP CONFIRMING FINDINGS (regen > scar, p < 0.05) ===")
sig = stats_df[(stats_df["p_value"] < 0.05) & (stats_df["log2fc_mean"] > 0)] \
      .sort_values("log2fc_mean", ascending=False)
print(sig.head(20).to_string(index=False))

print("\n=== TOP REFUTING FINDINGS (regen < scar, p < 0.05) ===")
refute = stats_df[(stats_df["p_value"] < 0.05) & (stats_df["log2fc_mean"] < 0)] \
      .sort_values("log2fc_mean")
print(refute.head(10).to_string(index=False))

print("\n=== DONE ===")
