"""
10 -- Compose 3 publication-ready Main Figures from panel-mode PNGs.

Each source PNG is expected to have NO internal suptitle and NO internal
a/b panel labels. This script adds:
  - one big "Main Figure N" title
  - one consistent set of external panel labels (a, b, c, d, ...)

Layout:
  MAIN FIGURE 1 -- Integration overview
    a. Fib scANVI Species UMAP
    b. Fib scANVI Leiden clusters
    c. KC scVI Species UMAP
    d. KC scVI Leiden clusters
    e. scIB benchmark bar chart
    f. Cluster × Species crosstab

  MAIN FIGURE 2 -- Fibroblast subtype identification
    a. Fib top-5 marker heatmap
    b. Fib canonical marker dotplot

  MAIN FIGURE 3 -- Cross-species hero gene validation
    a. Hero gene per-species dotplot
    b. Hero gene module-score UMAP + species boxplot
    c. Per-hero-gene expression box-plot (4 species)
"""
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from pathlib import Path

preview = Path("/Users/apple/Downloads/preview_output")
out_dir = preview / "composed_main_figures"
out_dir.mkdir(exist_ok=True, parents=True)

plt.rcParams.update({
    "font.family":       "Helvetica",
    "figure.facecolor":  "white",
    "savefig.facecolor": "white",
})


def load(name):
    p = preview / name
    if not p.exists():
        print(f"⚠️  missing: {p}")
        return None
    return mpimg.imread(p)


def add_panel_label(ax, label, fontsize=17):
    """Big bold a/b/c/... at the top-left of the axes."""
    ax.text(-0.01, 1.02, label, transform=ax.transAxes,
            fontsize=fontsize, fontweight="bold",
            va="bottom", ha="left")


# ==========================================================================
# MAIN FIGURE 1 — Integration overview
# ==========================================================================
print(">>> Composing Main Figure 1 — Integration overview ...")

fib_umap  = load("PREVIEW_CrossSpecies_Fib_UMAP_scANVI.png")
# prefer scANVI KC if available, fall back to scVI KC
kc_umap = load("PREVIEW_CrossSpecies_KC_UMAP_scANVI.png")
if kc_umap is None:
    kc_umap = load("PREVIEW_CrossSpecies_KC_UMAP.png")
scib_bar  = load("PREVIEW_scib_benchmark.png")
crosstab  = load("PREVIEW_cluster_species_crosstab.png")

fig = plt.figure(figsize=(15, 20))
gs  = fig.add_gridspec(4, 1, height_ratios=[1.15, 1.15, 0.75, 0.65],
                       hspace=0.08, left=0.04, right=0.98,
                       top=0.965, bottom=0.02)

ax1 = fig.add_subplot(gs[0, 0])
if fib_umap is not None: ax1.imshow(fib_umap)
ax1.axis("off")
add_panel_label(ax1, "a")
ax1.set_title("Fibroblast — scANVI integration",
              fontsize=12, fontweight="bold", pad=4, loc="left", x=0.02)

ax2 = fig.add_subplot(gs[1, 0])
if kc_umap is not None: ax2.imshow(kc_umap)
ax2.axis("off")
add_panel_label(ax2, "b")
ax2.set_title("Keratinocyte — scANVI integration",
              fontsize=12, fontweight="bold", pad=4, loc="left", x=0.02)

ax3 = fig.add_subplot(gs[2, 0])
if scib_bar is not None: ax3.imshow(scib_bar)
ax3.axis("off")
add_panel_label(ax3, "c")
ax3.set_title("scIB benchmark — scVI vs scANVI",
              fontsize=12, fontweight="bold", pad=4, loc="left", x=0.02)

ax4 = fig.add_subplot(gs[3, 0])
if crosstab is not None: ax4.imshow(crosstab)
ax4.axis("off")
add_panel_label(ax4, "d")
ax4.set_title("Fibroblast cluster × Species composition",
              fontsize=12, fontweight="bold", pad=4, loc="left", x=0.02)

fig.suptitle("Main Figure 1 · Cross-species integration overview",
             y=0.995, fontsize=18, fontweight="bold")
fig.savefig(out_dir / "MAIN_FIG_1_integration_overview.png",
            dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"    → {out_dir}/MAIN_FIG_1_integration_overview.png")

# ==========================================================================
# MAIN FIGURE 2 — Fibroblast subtype identification
# ==========================================================================
print("\n>>> Composing Main Figure 2 — Fib marker analysis ...")
fib_hm  = load("PREVIEW_Fib_marker_heatmap.png")
fib_dot = load("PREVIEW_Fib_marker_dotplot.png")

fig = plt.figure(figsize=(16, 12))
gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.15],
                      wspace=0.05, left=0.03, right=0.98,
                      top=0.94, bottom=0.03)

ax1 = fig.add_subplot(gs[0, 0])
if fib_hm is not None: ax1.imshow(fib_hm)
ax1.axis("off")
add_panel_label(ax1, "a")
ax1.set_title("Top 5 unbiased markers per cluster (Wilcoxon)",
              fontsize=12, fontweight="bold", pad=4, loc="left", x=0.02)

ax2 = fig.add_subplot(gs[0, 1])
if fib_dot is not None: ax2.imshow(fib_dot)
ax2.axis("off")
add_panel_label(ax2, "b")
ax2.set_title("Canonical fibroblast subtype markers",
              fontsize=12, fontweight="bold", pad=4, loc="left", x=0.02)

fig.suptitle("Main Figure 2 · Fibroblast subtype identification",
             y=0.99, fontsize=18, fontweight="bold")
fig.savefig(out_dir / "MAIN_FIG_2_fib_markers.png",
            dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"    → {out_dir}/MAIN_FIG_2_fib_markers.png")

# ==========================================================================
# MAIN FIGURE 3 — Hero gene validation
# ==========================================================================
print("\n>>> Composing Main Figure 3 — Hero gene validation ...")
hero_dot   = load("PREVIEW_hero_dotplot_species.png")
hero_umap  = load("PREVIEW_hero_module_score_umap.png")
hero_box   = load("PREVIEW_hero_boxplot_by_species.png")

fig = plt.figure(figsize=(16, 15))
gs = fig.add_gridspec(3, 1, height_ratios=[1.2, 0.9, 0.85],
                       hspace=0.12, left=0.03, right=0.98,
                       top=0.955, bottom=0.03)

ax1 = fig.add_subplot(gs[0, 0])
if hero_dot is not None: ax1.imshow(hero_dot)
ax1.axis("off")
add_panel_label(ax1, "a")
ax1.set_title("Hero gene expression per cluster, split by species",
              fontsize=12, fontweight="bold", pad=4, loc="left", x=0.02)

ax2 = fig.add_subplot(gs[1, 0])
if hero_umap is not None: ax2.imshow(hero_umap)
ax2.axis("off")
add_panel_label(ax2, "b")
ax2.set_title("Hero-gene module score across UMAP and species",
              fontsize=12, fontweight="bold", pad=4, loc="left", x=0.02)

ax3 = fig.add_subplot(gs[2, 0])
if hero_box is not None: ax3.imshow(hero_box)
ax3.axis("off")
add_panel_label(ax3, "c")
ax3.set_title("Per-gene expression across 4 species",
              fontsize=12, fontweight="bold", pad=4, loc="left", x=0.02)

fig.suptitle("Main Figure 3 · Cross-species hero gene validation",
             y=0.99, fontsize=18, fontweight="bold")
fig.savefig(out_dir / "MAIN_FIG_3_hero_validation.png",
            dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"    → {out_dir}/MAIN_FIG_3_hero_validation.png")

print("\n=== DONE ===")
print(f"\nAll 3 main figures written to:\n  {out_dir}/")
