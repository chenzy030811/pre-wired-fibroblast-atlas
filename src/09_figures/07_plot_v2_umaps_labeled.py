"""
Plot v2 UMAPs with BIOLOGICAL LABELS instead of just cluster numbers.

Labels based on marker analysis (05d fib, 05e kc).
Output: annotated versions of Fib scANVI and KC scVI UMAPs.
"""
import anndata as ad
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

base    = Path("/Users/apple/Downloads/研究生毕业论文/Cross_Species")
out_dir = Path("/Users/apple/Downloads/preview_output")

# ================ biological labels ================
# Fibroblast scANVI (14 clusters, after reclustering)
FIB_LABELS = {
    "0":  "0 · Acomys-Fib",
    "1":  "1 · Progenitor (GREM1+)",
    "2":  "2 · Mus-Fib",
    "3":  "3 · Papillary",
    "4":  "4 · Universal-Activated",
    "5":  "5 · WNT/adipogenic",
    "6":  "6 · Reticular",
    "7":  "7 · Human-Fib",
    "8":  "8 · Human-Fib",
    "9":  "9 · Pre-wired★",
    "10": "10 · CHRDL1-high",
    "11": "11 · WNT-niche (RSPO3+)",
    "12": "12 · Human-Fib",
    "13": "13 · Myofibroblast",
}

# Keratinocyte scVI (12 clusters after reclustering)
KC_LABELS = {
    "0":  "0 · HFSC (LGR5+)",
    "1":  "1 · HF-related",
    "2":  "2 · Outer root sheath",
    "3":  "3 · Basal",
    "4":  "4 · Suprabasal (early)",
    "5":  "5 · Basal-like",
    "6":  "6 · Cycling",
    "7":  "7 · Wound-response",
    "8":  "8 · HF outer",
    "9":  "9 · Sebocyte",
    "10": "10 · Suprabasal (mature)",
    "11": "11 · Wound / HF",
}

SPECIES_COLORS = {
    "Reindeer_Antler": "#E28A7A",
    "Reindeer_Back":   "#8A5A55",
    "Human":           "#8AB2D0",
    "Mus":             "#B5C275",
    "Acomys":          "#D189A5",
}
LEIDEN_COLORS = [
    "#E28A7A", "#D4B14A", "#E4A08C", "#7BC5AF",
    "#B5C275", "#A98CBE", "#8AB2D0", "#6B8DBA",
    "#D189A5", "#B387A2", "#9BA0AE", "#C29A6E",
    "#8FD1BB", "#B989A2",
]

plt.rcParams.update({
    "font.family":"Helvetica", "font.size":11,
    "axes.titlesize":14, "axes.titleweight":"bold",
    "figure.facecolor":"white", "savefig.facecolor":"white",
})

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

def plot_labeled(in_path, leiden_key, labels_dict, title, out_path):
    print(f">>> {out_path.name} ...")
    adata = ad.read_h5ad(in_path)
    umap = adata.obsm["X_umap"]
    species = adata.obs["species"].values.astype(str)
    leiden  = adata.obs[leiden_key].values.astype(str)

    fig, axes = plt.subplots(1, 2, figsize=(15, 6.5), dpi=200)
    xr = umap[:,0].min()-1, umap[:,0].max()+1
    yr = umap[:,1].min()-1, umap[:,1].max()+1

    # panel a: species
    ax = axes[0]
    for sp in ["Reindeer_Antler","Reindeer_Back","Human","Mus","Acomys"]:
        m = species == sp
        if m.sum() == 0: continue
        ax.scatter(umap[m,0], umap[m,1], s=2.5, alpha=0.6,
                   c=SPECIES_COLORS[sp], label=f"{sp} (n={m.sum():,})",
                   edgecolors="none", rasterized=True)
    ax.set_xlim(xr); ax.set_ylim(yr)
    strip_axes(ax); corner_axes(ax, xr, yr)
    ax.set_title("Species", pad=8)
    _handles, _labels = ax.get_legend_handles_labels()

    # panel b: leiden with biological labels
    ax = axes[1]
    uniq = sorted(np.unique(leiden), key=int)
    for i, c in enumerate(uniq):
        m = leiden == c
        ax.scatter(umap[m,0], umap[m,1], s=2.5, alpha=0.6,
                   c=LEIDEN_COLORS[i % len(LEIDEN_COLORS)],
                   edgecolors="none", rasterized=True)
        if m.sum() > 0:
            cx = np.median(umap[m,0])
            cy = np.median(umap[m,1])
            label = labels_dict.get(c, c)
            # ★ star for key hero-gene clusters
            fontweight = "bold" if "★" in label else "medium"
            fontsize   = 8.5 if "★" in label else 7.5
            color      = "#C63C4B" if "★" in label else "#222222"
            ax.text(cx, cy, label, fontsize=fontsize, fontweight=fontweight,
                    ha="center", va="center", color=color,
                    family="Helvetica", zorder=10,
                    bbox=dict(boxstyle="round,pad=0.15",
                              facecolor="white", edgecolor="none", alpha=0.85))
    ax.set_xlim(xr); ax.set_ylim(yr)
    strip_axes(ax); corner_axes(ax, xr, yr)
    ax.set_title("Fibroblast subtypes"
                 if "Fib" in title else "Keratinocyte subtypes", pad=8)

    leg = fig.legend(_handles, _labels,
                     loc="lower center", bbox_to_anchor=(0.5, -0.02),
                     ncol=5, frameon=False, fontsize=9,
                     markerscale=4, handletextpad=0.4, columnspacing=1.2)
    for lh in leg.legend_handles: lh.set_alpha(1.0)

    plt.suptitle(title, y=0.99, fontsize=13, fontweight="bold")
    plt.tight_layout(rect=[0, 0.05, 1, 0.96])
    plt.savefig(out_path, dpi=600, bbox_inches="tight")
    plt.close()
    print(f"    → {out_path}")

plot_labeled(
    base/"scvi_out/Fibroblast_v2_scanvi/integrated.h5ad",
    "leiden_scanvi", FIB_LABELS,
    "Cross-species Fibroblast integration — scANVI · with subtype labels",
    out_dir/"PREVIEW_v2_Fib_UMAP_labeled.png"
)
plot_labeled(
    base/"scvi_out/Keratinocyte_v2_scvi/integrated.h5ad",
    "leiden_scvi", KC_LABELS,
    "Cross-species Keratinocyte integration — scVI · with subtype labels",
    out_dir/"PREVIEW_v2_KC_UMAP_labeled.png"
)
print("=== DONE ===")
