"""
09 -- Summary comparison figure: scVI vs scANVI integration.

Layout (single figure):
  Row 1: scVI    — [Species UMAP] | [Leiden UMAP]
  Row 2: scANVI  — [Species UMAP] | [Leiden UMAP]
  Row 3: scIB benchmark bar chart with all metrics

One publication-quality figure that captures the full comparison
+ the scIB verdict.
"""
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# =========== paths ===========
base    = Path("/Users/apple/Downloads/研究生毕业论文/Cross_Species")
scvi_h5  = base / "scvi_out/Fibroblast_scvi/integrated.h5ad"
scanvi_h5= base / "scvi_out/Fibroblast_scanvi/integrated.h5ad"
scib_csv = base / "results/scib_benchmark.csv"
preview  = Path("/Users/apple/Downloads/preview_output")

# =========== palettes ===========
SPECIES_COLORS = {
    "Reindeer": "#C66055",
    "Human":    "#5090C6",
    "Mus":      "#7C8E2C",
    "Acomys":   "#B83668",
}
LEIDEN_COLORS = [
    "#C66055", "#A07C20", "#AD5818", "#2D9686",
    "#7C8E2C", "#9560AE", "#5090C6", "#1F4F94",
    "#B83668", "#85316C", "#4A4E69", "#6D4C1F",
    "#7BC5AF", "#A96F8C", "#B89364", "#6082B6",
]

plt.rcParams.update({
    "font.family":     "Helvetica",
    "font.size":       10,
    "axes.linewidth":  0.6,
    "figure.facecolor":"white",
    "savefig.facecolor":"white",
})

# =========== load ===========
print(">>> Loading scVI + scANVI + scIB ...")
a_sv = ad.read_h5ad(scvi_h5);   a_sv.obs_names_make_unique()
a_sa = ad.read_h5ad(scanvi_h5); a_sa.obs_names_make_unique()
scib = pd.read_csv(scib_csv, index_col=0)
print(f"    scVI:   {a_sv.n_obs} cells, {a_sv.obs['leiden_scvi'].nunique()} clusters")
print(f"    scANVI: {a_sa.n_obs} cells, {a_sa.obs['leiden_scanvi'].nunique()} clusters")

# =========== helper ===========
def corner_axes(ax, xr, yr, frac=0.14, gap=0.02):
    xs = xr[0]+(xr[1]-xr[0])*gap; ys = yr[0]+(yr[1]-yr[0])*gap
    xl = xs+(xr[1]-xr[0])*frac;   yl = ys+(yr[1]-yr[0])*frac
    ax.annotate("", xy=(xl,ys), xytext=(xs,ys),
                arrowprops=dict(arrowstyle="->", color="black", lw=0.8))
    ax.annotate("", xy=(xs,yl), xytext=(xs,ys),
                arrowprops=dict(arrowstyle="->", color="black", lw=0.8))
    ax.text((xs+xl)/2, ys-(yr[1]-yr[0])*0.03, "UMAP 1",
            ha="center", va="top", fontsize=7)
    ax.text(xs-(xr[1]-xr[0])*0.02, (ys+yl)/2, "UMAP 2",
            ha="right", va="center", fontsize=7, rotation=90)

def strip_axes(ax):
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_visible(False)

def plot_species(ax, umap, species):
    xr = umap[:,0].min()-1, umap[:,0].max()+1
    yr = umap[:,1].min()-1, umap[:,1].max()+1
    for sp, c in SPECIES_COLORS.items():
        m = species == sp
        ax.scatter(umap[m,0], umap[m,1], s=1.2, alpha=0.55, c=c,
                   label=f"{sp} (n={m.sum():,})",
                   edgecolors="none", rasterized=True)
    ax.set_xlim(xr); ax.set_ylim(yr)
    strip_axes(ax); corner_axes(ax, xr, yr)
    leg = ax.legend(loc="upper right", frameon=False, fontsize=7,
                    markerscale=4, handletextpad=0.2, borderaxespad=0.2)
    for lh in leg.legend_handles: lh.set_alpha(1.0)

def plot_leiden(ax, umap, leiden):
    xr = umap[:,0].min()-1, umap[:,0].max()+1
    yr = umap[:,1].min()-1, umap[:,1].max()+1
    uniq = sorted(np.unique(leiden), key=int)
    for i, c in enumerate(uniq):
        m = leiden == c
        ax.scatter(umap[m,0], umap[m,1], s=1.2, alpha=0.55,
                   c=LEIDEN_COLORS[i % len(LEIDEN_COLORS)],
                   edgecolors="none", rasterized=True)
        if m.sum() > 0:
            cx = np.median(umap[m,0]); cy = np.median(umap[m,1])
            ax.text(cx, cy, c, fontsize=7, fontweight="bold",
                    ha="center", va="center", color="black",
                    bbox=dict(boxstyle="circle,pad=0.18",
                              facecolor="white", edgecolor="black", linewidth=0.4))
    ax.set_xlim(xr); ax.set_ylim(yr)
    strip_axes(ax); corner_axes(ax, xr, yr)

# =========== figure layout ===========
fig = plt.figure(figsize=(13, 14))
gs  = fig.add_gridspec(3, 2, height_ratios=[1, 1, 0.9],
                       hspace=0.20, wspace=0.08,
                       left=0.05, right=0.98, top=0.94, bottom=0.05)

# ----- Row 1: scVI -----
ax_a = fig.add_subplot(gs[0, 0])
plot_species(ax_a, a_sv.obsm["X_umap"], a_sv.obs["species"].values.astype(str))
ax_a.set_title("scVI — Species", fontsize=12, fontweight="bold", pad=6)

ax_b = fig.add_subplot(gs[0, 1])
plot_leiden(ax_b, a_sv.obsm["X_umap"], a_sv.obs["leiden_scvi"].values.astype(str))
ax_b.set_title(f"scVI — Leiden clusters (k={a_sv.obs['leiden_scvi'].nunique()})",
               fontsize=12, fontweight="bold", pad=6)

# ----- Row 2: scANVI -----
ax_c = fig.add_subplot(gs[1, 0])
plot_species(ax_c, a_sa.obsm["X_umap"], a_sa.obs["species"].values.astype(str))
ax_c.set_title("scANVI — Species", fontsize=12, fontweight="bold", pad=6)

ax_d = fig.add_subplot(gs[1, 1])
plot_leiden(ax_d, a_sa.obsm["X_umap"], a_sa.obs["leiden_scanvi"].values.astype(str))
ax_d.set_title(f"scANVI — Leiden clusters (k={a_sa.obs['leiden_scanvi'].nunique()})",
               fontsize=12, fontweight="bold", pad=6)

# ----- Row 3: scIB benchmark bar chart -----
ax_e = fig.add_subplot(gs[2, :])
# select numeric metric rows only (drop text/summary)
metric_rows = [i for i in scib.index if i not in {"Metric Type"}]
plot_df = scib.loc[metric_rows].copy()
# separate main aggregates from individual metrics
summary_cols = ["Bio conservation", "Batch correction", "Total"]
main_metrics = [c for c in plot_df.columns if c not in summary_cols]

# get scVI vs scANVI values
scvi_vals = plot_df.loc[[i for i in ["X_scVI"] if i in plot_df.index][0], main_metrics].astype(float).values
scanvi_vals = plot_df.loc[[i for i in ["X_scANVI"] if i in plot_df.index][0], main_metrics].astype(float).values

xpos = np.arange(len(main_metrics))
w = 0.36
ax_e.bar(xpos-w/2, scvi_vals, w, label="scVI",
         color="#5090C6", edgecolor="black", linewidth=0.4)
ax_e.bar(xpos+w/2, scanvi_vals, w, label="scANVI",
         color="#C66055", edgecolor="black", linewidth=0.4)

# annotate values above bars
for i, v in enumerate(scvi_vals):
    ax_e.text(i-w/2, v+0.015, f"{v:.2f}", ha="center", va="bottom", fontsize=8)
for i, v in enumerate(scanvi_vals):
    ax_e.text(i+w/2, v+0.015, f"{v:.2f}", ha="center", va="bottom", fontsize=8)

ax_e.set_xticks(xpos)
ax_e.set_xticklabels(main_metrics, rotation=25, ha="right", fontsize=9)
ax_e.set_ylabel("Score (higher = better)", fontsize=10)
ax_e.set_ylim(0, 1.08)
ax_e.set_title("scIB benchmark — quantitative comparison "
               f"(Total: scVI={float(scib.loc['X_scVI', 'Total']):.3f}  "
               f"vs  scANVI={float(scib.loc['X_scANVI', 'Total']):.3f}  →  scANVI wins)",
               fontsize=12, fontweight="bold", pad=8)
ax_e.legend(frameon=False, fontsize=10, loc="lower right")
for s in ["top", "right"]: ax_e.spines[s].set_visible(False)
ax_e.grid(axis="y", linestyle=":", linewidth=0.4, color="#D0D0D0")
ax_e.set_axisbelow(True)

# ----- panel labels -----
for ax_i, label in zip([ax_a, ax_b, ax_c, ax_d, ax_e],
                        ["a", "b", "c", "d", "e"]):
    ax_i.text(-0.02, 1.05, label, transform=ax_i.transAxes,
              fontsize=14, fontweight="bold", va="bottom", ha="right")

fig.suptitle("Cross-species Fibroblast integration — scVI vs scANVI comparison",
             y=0.985, fontsize=14, fontweight="bold")

out = preview / "PREVIEW_scvi_vs_scanvi_summary.png"
plt.savefig(out, dpi=400, bbox_inches="tight")
plt.close()
print(f"\nWrote {out}")
