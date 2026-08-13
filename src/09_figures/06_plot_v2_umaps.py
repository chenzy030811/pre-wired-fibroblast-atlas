"""
Plot 4 UMAPs for v2 integrated data (5 groups):
  - Fibroblast scVI    (loser reference)
  - Fibroblast scANVI  (WINNER — used downstream)
  - Keratinocyte scVI  (WINNER — used downstream)
  - Keratinocyte scANVI (loser reference)

Unified pastel style, bottom legend, plain black cluster numbers.
"""
import anndata as ad
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# =========== paths ===========
base = Path("/Users/apple/Downloads/研究生毕业论文/Cross_Species")
out_dir = Path("/Users/apple/Downloads/preview_output")

# =========== unified palette (pastel, matches all main figures) ===========
SPECIES_COLORS = {
    "Reindeer_Antler": "#D4776A",   # deeper coral — regen
    "Reindeer_Back":   "#8A5A55",   # brownish coral — scar
    "Human":           "#8AB2D0",
    "Mus":             "#B5C275",
    "Acomys":           "#D189A5",
}
LEIDEN_COLORS = [
    "#E28A7A", "#D4B14A", "#E4A08C", "#7BC5AF",
    "#B5C275", "#A98CBE", "#8AB2D0", "#6B8DBA",
    "#D189A5", "#B387A2", "#9BA0AE", "#C29A6E",
    "#8FD1BB", "#B989A2", "#C7A67E", "#7C97C6",
    "#A6D8C7", "#D19FB4", "#B9A56F", "#95A5C6",
]

plt.rcParams.update({
    "font.family":       "Helvetica",
    "font.size":         11,
    "axes.titlesize":    14,
    "axes.titleweight":  "bold",
    "figure.facecolor":  "white",
    "savefig.facecolor": "white",
})

# =========== helpers ===========
def corner_axes(ax, xr, yr, frac=0.14, gap=0.02):
    xs = xr[0]+(xr[1]-xr[0])*gap; ys = yr[0]+(yr[1]-yr[0])*gap
    xl = xs+(xr[1]-xr[0])*frac;   yl = ys+(yr[1]-yr[0])*frac
    ax.annotate("", xy=(xl,ys), xytext=(xs,ys),
                arrowprops=dict(arrowstyle="->", color="black", lw=1.0))
    ax.annotate("", xy=(xs,yl), xytext=(xs,ys),
                arrowprops=dict(arrowstyle="->", color="black", lw=1.0))
    ax.text((xs+xl)/2, ys-(yr[1]-yr[0])*0.03, "UMAP 1",
            ha="center", va="top", fontsize=9)
    ax.text(xs-(xr[1]-xr[0])*0.02, (ys+yl)/2, "UMAP 2",
            ha="right", va="center", fontsize=9, rotation=90)

def strip_axes(ax):
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_visible(False)

def plot_one(in_path, leiden_key, method_label, out_path):
    print(f">>> {out_path.name} ...")
    adata = ad.read_h5ad(in_path)
    umap = adata.obsm["X_umap"]
    species = adata.obs["species"].values.astype(str)
    leiden  = adata.obs[leiden_key].values.astype(str)
    n_clust = len(np.unique(leiden))
    print(f"    {adata.n_obs} cells · {n_clust} clusters ({method_label})")

    fig, axes = plt.subplots(1, 2, figsize=(13, 6.2), dpi=200)
    xr = umap[:,0].min()-1, umap[:,0].max()+1
    yr = umap[:,1].min()-1, umap[:,1].max()+1

    # ---- Panel a: species (5 groups) ----
    ax = axes[0]
    for sp in ["Reindeer_Antler", "Reindeer_Back", "Human", "Mus", "Acomys"]:
        m = species == sp
        if m.sum() == 0: continue
        ax.scatter(umap[m,0], umap[m,1],
                   s=2.5, alpha=0.6, c=SPECIES_COLORS[sp],
                   label=f"{sp} (n={m.sum():,})",
                   edgecolors="none", linewidths=0, rasterized=True)
    ax.set_xlim(xr); ax.set_ylim(yr)
    strip_axes(ax); corner_axes(ax, xr, yr)
    ax.set_title("Species", pad=8)
    _handles, _labels = ax.get_legend_handles_labels()

    # ---- Panel b: leiden clusters ----
    ax = axes[1]
    uniq = sorted(np.unique(leiden), key=int)
    for i, c in enumerate(uniq):
        m = leiden == c
        ax.scatter(umap[m,0], umap[m,1],
                   s=2.5, alpha=0.6,
                   c=LEIDEN_COLORS[i % len(LEIDEN_COLORS)],
                   edgecolors="none", linewidths=0, rasterized=True)
        if m.sum() > 0:
            cx = np.median(umap[m,0])
            cy = np.median(umap[m,1])
            ax.text(cx, cy, c, fontsize=9, fontweight="medium",
                    ha="center", va="center", color="#222222",
                    family="Helvetica", zorder=10)
    ax.set_xlim(xr); ax.set_ylim(yr)
    strip_axes(ax); corner_axes(ax, xr, yr)
    ax.set_title(f"Leiden clusters ({method_label})", pad=8)

    # figure-level species legend at bottom (5 items)
    leg = fig.legend(_handles, _labels,
                     loc="lower center", bbox_to_anchor=(0.5, -0.02),
                     ncol=5, frameon=False, fontsize=9,
                     markerscale=4, handletextpad=0.4, columnspacing=1.2)
    for lh in leg.legend_handles: lh.set_alpha(1.0)

    plt.tight_layout(rect=[0, 0.05, 1, 0.98])
    plt.savefig(out_path, dpi=600, bbox_inches="tight")
    plt.close()
    print(f"    → {out_path}")

# =========== run all 4 ===========
plot_one(base/"scvi_out/Fibroblast_v2_scvi/integrated.h5ad",     "leiden_scvi",
         "scVI",   out_dir/"PREVIEW_v2_Fib_scVI_UMAP.png")
plot_one(base/"scvi_out/Fibroblast_v2_scanvi/integrated.h5ad",   "leiden_scanvi",
         "scANVI", out_dir/"PREVIEW_v2_Fib_scANVI_UMAP.png")
plot_one(base/"scvi_out/Keratinocyte_v2_scvi/integrated.h5ad",   "leiden_scvi",
         "scVI",   out_dir/"PREVIEW_v2_KC_scVI_UMAP.png")
plot_one(base/"scvi_out/Keratinocyte_v2_scanvi/integrated.h5ad", "leiden_scanvi",
         "scANVI", out_dir/"PREVIEW_v2_KC_scANVI_UMAP.png")

print("\n=== DONE ===")
