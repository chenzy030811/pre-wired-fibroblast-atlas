"""
Publication-grade UMAP for cross-species Fibroblast integration.
Cloud-style points + Cell-journal palette + corner UMAP1/UMAP2 arrow.

Uses REAL field names verified from integrated.h5ad:
  .obsm['X_umap']         2D UMAP coords
  .obs['species']         4 categories: Reindeer / Human / Mus / Acomys
  .obs['leiden_scvi']     12 clusters (0..11)
"""
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path

# =========== paths ===========
in_path  = Path("/Users/apple/Downloads/研究生毕业论文/Cross_Species/"
                "scvi_out/Fibroblast_scvi/integrated.h5ad")
out_dir  = Path("/Users/apple/Downloads/preview_output")
out_dir.mkdir(exist_ok=True, parents=True)

# =========== palettes — Reynolds/Ferreira pastel, slightly saturated ===========
SPECIES_COLORS = {
    "Reindeer": "#D4776A",   # coral (a touch deeper)
    "Human":    "#6B95C4",   # sky blue (a touch deeper)
    "Mus":      "#9EAE60",   # olive (a touch deeper)
    "Acomys":   "#C0779A",   # rose (a touch deeper)
}
LEIDEN_COLORS = [
    "#D4776A", "#C09A38", "#D68F79", "#5DB598",
    "#9EAE60", "#9377AE", "#6B95C4", "#4F70A6",
    "#C0779A", "#95607D", "#83889A", "#B0885B",
]

# =========== load ===========
print(f"Loading {in_path.name} ...")
adata = ad.read_h5ad(in_path)
umap = adata.obsm["X_umap"]
species  = adata.obs["species"].values.astype(str)
leiden   = adata.obs["leiden_scvi"].values.astype(str)
print(f"  {adata.n_obs} cells")

# =========== typography ===========
plt.rcParams.update({
    "font.family":     "Helvetica",
    "font.size":       11,
    "axes.titlesize":  14,
    "axes.titleweight":"bold",
    "figure.facecolor":"white",
    "savefig.facecolor":"white",
})

# =========== helper: corner UMAP1/UMAP2 arrow  (like reindeer figure) ===========
def corner_axes(ax, xr, yr, frac=0.14, gap=0.02):
    """Draw a small L-shaped arrow indicator at bottom-left."""
    xs = xr[0] + (xr[1]-xr[0]) * gap
    ys = yr[0] + (yr[1]-yr[0]) * gap
    xl = xs + (xr[1]-xr[0]) * frac
    yl = ys + (yr[1]-yr[0]) * frac
    # horizontal arrow -> UMAP1
    ax.annotate("", xy=(xl, ys), xytext=(xs, ys),
                arrowprops=dict(arrowstyle="->", color="black", lw=1.0))
    # vertical arrow -> UMAP2
    ax.annotate("", xy=(xs, yl), xytext=(xs, ys),
                arrowprops=dict(arrowstyle="->", color="black", lw=1.0))
    ax.text((xs+xl)/2, ys - (yr[1]-yr[0])*0.03,
            "UMAP 1", ha="center", va="top", fontsize=9)
    ax.text(xs - (xr[1]-xr[0])*0.02, (ys+yl)/2,
            "UMAP 2", ha="right", va="center", fontsize=9, rotation=90)


def strip_axes(ax):
    """Remove all axis lines / ticks — corner arrow shows orientation."""
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


# =========== figure ===========
fig, axes = plt.subplots(1, 2, figsize=(13, 6.2), dpi=200)

xr = umap[:,0].min()-1, umap[:,0].max()+1
yr = umap[:,1].min()-1, umap[:,1].max()+1

# ---------- Panel 1: species ----------
ax = axes[0]
for sp, color in SPECIES_COLORS.items():
    mask = species == sp
    ax.scatter(umap[mask,0], umap[mask,1],
               s=2.5, alpha=0.6, c=color, label=f"{sp} (n={mask.sum():,})",
               edgecolors="none", linewidths=0, rasterized=True)
ax.set_xlim(xr); ax.set_ylim(yr)
strip_axes(ax)
corner_axes(ax, xr, yr)
ax.set_title("Species", pad=8)
leg = ax.legend(loc="upper right", frameon=False, fontsize=9,
                markerscale=3, handletextpad=0.2, borderaxespad=0.2)
for lh in leg.legend_handles:
    lh.set_alpha(1.0)

# ---------- Panel 2: leiden clusters ----------
ax = axes[1]
for i in range(12):
    mask = leiden == str(i)
    ax.scatter(umap[mask,0], umap[mask,1],
               s=2.5, alpha=0.6, c=LEIDEN_COLORS[i],
               edgecolors="none", linewidths=0, rasterized=True)
    # plain black number label at cluster centre — no circle, no halo
    if mask.sum() > 0:
        cx = np.median(umap[mask,0])
        cy = np.median(umap[mask,1])
        # Reynolds style: plain black text, no box, just clean label
        ax.text(cx, cy, str(i), fontsize=9, fontweight="medium",
                ha="center", va="center", color="#222222",
                family="Helvetica")
ax.set_xlim(xr); ax.set_ylim(yr)
strip_axes(ax)
corner_axes(ax, xr, yr)
ax.set_title("Leiden clusters (scVI)", pad=8)

plt.suptitle(
    f"Cross-species Fibroblast integration — {adata.n_obs:,} cells, "
    f"scVI (200 epochs)",
    y=0.99, fontsize=13, fontweight="bold"
)

plt.tight_layout(rect=[0,0,1,0.96])

out_path = out_dir / "PREVIEW_CrossSpecies_Fib_UMAP.png"
plt.savefig(out_path, dpi=600, bbox_inches="tight")
print(f"\nWrote {out_path}")
