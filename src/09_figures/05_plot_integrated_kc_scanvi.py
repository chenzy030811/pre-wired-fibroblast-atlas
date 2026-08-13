"""
Publication-grade UMAP for cross-species Keratinocyte — scANVI integration.
Same style as Fib scANVI (Reynolds pastel, cloud, bottom legend, no chrome).

Fields read from scanvi integrated.h5ad:
  .obsm['X_umap']         2D UMAP from scANVI latent
  .obs['species']         4 species
  .obs['leiden_scanvi']   new cluster assignments after scANVI
"""
import anndata as ad
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# =========== paths ===========
in_path  = Path("/Users/apple/Downloads/研究生毕业论文/Cross_Species/"
                "scvi_out/Keratinocyte_scanvi/integrated.h5ad")
out_dir  = Path("/Users/apple/Downloads/preview_output")
out_dir.mkdir(exist_ok=True, parents=True)

# =========== palette — unified pastel across all main figures ==
SPECIES_COLORS = {
    "Reindeer": "#E28A7A",
    "Human":    "#8AB2D0",
    "Mus":      "#B5C275",
    "Acomys":   "#D189A5",
}
LEIDEN_COLORS = [
    "#E28A7A", "#D4B14A", "#E4A08C", "#7BC5AF",
    "#B5C275", "#A98CBE", "#8AB2D0", "#6B8DBA",
    "#D189A5", "#B387A2", "#9BA0AE", "#C29A6E",
    "#8FD1BB", "#B989A2", "#C7A67E", "#7C97C6",
    "#A6D8C7", "#D19FB4", "#B9A56F",
]

# =========== load ===========
print(f"Loading {in_path.name} ...")
adata = ad.read_h5ad(in_path)
umap = adata.obsm["X_umap"]
species = adata.obs["species"].values.astype(str)
leiden  = adata.obs["leiden_scanvi"].values.astype(str)
n_clust = len(np.unique(leiden))
print(f"  {adata.n_obs} cells · {n_clust} clusters (scANVI)")

# =========== typography ===========
plt.rcParams.update({
    "font.family":     "Helvetica",
    "font.size":       11,
    "axes.titlesize":  14,
    "axes.titleweight":"bold",
    "figure.facecolor":"white",
    "savefig.facecolor":"white",
})

def corner_axes(ax, xr, yr, frac=0.14, gap=0.02):
    xs = xr[0] + (xr[1]-xr[0]) * gap
    ys = yr[0] + (yr[1]-yr[0]) * gap
    xl = xs + (xr[1]-xr[0]) * frac
    yl = ys + (yr[1]-yr[0]) * frac
    ax.annotate("", xy=(xl, ys), xytext=(xs, ys),
                arrowprops=dict(arrowstyle="->", color="black", lw=1.0))
    ax.annotate("", xy=(xs, yl), xytext=(xs, ys),
                arrowprops=dict(arrowstyle="->", color="black", lw=1.0))
    ax.text((xs+xl)/2, ys - (yr[1]-yr[0])*0.03,
            "UMAP 1", ha="center", va="top", fontsize=9)
    ax.text(xs - (xr[1]-xr[0])*0.02, (ys+yl)/2,
            "UMAP 2", ha="right", va="center", fontsize=9, rotation=90)

def strip_axes(ax):
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_visible(False)

# =========== figure ===========
fig, axes = plt.subplots(1, 2, figsize=(13, 6.2), dpi=200)

xr = umap[:,0].min()-1, umap[:,0].max()+1
yr = umap[:,1].min()-1, umap[:,1].max()+1

# ---- panel: species ----
ax = axes[0]
for sp, color in SPECIES_COLORS.items():
    mask = species == sp
    ax.scatter(umap[mask,0], umap[mask,1],
               s=2.5, alpha=0.6, c=color, label=f"{sp} (n={mask.sum():,})",
               edgecolors="none", linewidths=0, rasterized=True)
ax.set_xlim(xr); ax.set_ylim(yr)
strip_axes(ax); corner_axes(ax, xr, yr)
ax.set_title("Species", pad=8)
_handles, _labels = ax.get_legend_handles_labels()

# ---- panel: leiden clusters ----
ax = axes[1]
uniq_clust = sorted(np.unique(leiden), key=int)
for i, c in enumerate(uniq_clust):
    mask = leiden == c
    ax.scatter(umap[mask,0], umap[mask,1],
               s=2.5, alpha=0.6, c=LEIDEN_COLORS[i % len(LEIDEN_COLORS)],
               edgecolors="none", linewidths=0, rasterized=True)
    if mask.sum() > 0:
        cx = np.median(umap[mask,0])
        cy = np.median(umap[mask,1])
        ax.text(cx, cy, c, fontsize=9, fontweight="medium",
                ha="center", va="center", color="#222222",
                family="Helvetica", zorder=10)
ax.set_xlim(xr); ax.set_ylim(yr)
strip_axes(ax); corner_axes(ax, xr, yr)
ax.set_title("Leiden clusters (scANVI)", pad=8)

# standalone mode: with title, panel labels a/b, bottom legend
for ax_i, label in zip(axes, ["a", "b"]):
    ax_i.text(-0.02, 1.02, label, transform=ax_i.transAxes,
              fontsize=15, fontweight="bold", va="bottom", ha="right")

plt.suptitle(
    f"Cross-species Keratinocyte integration — scANVI (n = {adata.n_obs:,})",
    y=0.99, fontsize=13, fontweight="bold"
)

leg = fig.legend(_handles, _labels,
                 loc="lower center", bbox_to_anchor=(0.5, -0.02),
                 ncol=4, frameon=False, fontsize=10,
                 markerscale=4, handletextpad=0.4, columnspacing=1.8)
for lh in leg.legend_handles: lh.set_alpha(1.0)

plt.tight_layout(rect=[0, 0.04, 1, 0.96])
out_path = out_dir / "PREVIEW_CrossSpecies_KC_UMAP_scANVI.png"
plt.savefig(out_path, dpi=600, bbox_inches="tight")
print(f"\nWrote {out_path}")
