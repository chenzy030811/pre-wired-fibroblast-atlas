"""
Create a schematic 'Fig 2 · hero genes' concept diagram.
Shows 7 hero genes grouped by their signalling role
(WNT-active / BMP-antagonised / immunomodulatory / basement-membrane).
"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch

# palette
NAVY  = "#1B2A4E"
CORAL = "#C66055"
BLUE  = "#5090C6"
OLIVE = "#7C8E2C"
ROSE  = "#B83668"
CREAM = "#F5F1EB"
DARK  = "#222222"
MUTED = "#7A7A7A"

fig, ax = plt.subplots(figsize=(10, 5.5), dpi=200)
ax.set_xlim(0, 10); ax.set_ylim(0, 5.5); ax.axis("off")

# central niche label
ax.text(5, 5.0, "Pre-wired regenerative fibroblast niche",
        ha="center", va="top", fontsize=14, fontweight="bold",
        color=NAVY, family="Helvetica")
ax.text(5, 4.6, "constitutively expressed in Acomys dermal fibroblasts",
        ha="center", va="top", fontsize=10, style="italic", color=MUTED,
        family="Helvetica")

# 4 grouped boxes for functional categories
groups = [
    ("WNT-activating",       ["WNT2", "RSPO3"],       CORAL,  0.5, 2.3),
    ("BMP-antagonising",     ["GREM1", "SULF1", "CHRDL1"], BLUE, 3.3, 2.3),
    ("Inflammation buffer",  ["IL13RA2"],             OLIVE,  6.6, 2.3),
    ("Follicular ECM",       ["COL13A1"],             ROSE,   8.3, 2.3),
]

box_w = 1.5; box_h = 1.7
for title, genes, color, x, y in groups:
    if title == "BMP-antagonising":
        w_here = 2.9
    elif title == "WNT-activating":
        w_here = 2.5
    else:
        w_here = box_w
    r = FancyBboxPatch((x, y), w_here, box_h,
                       boxstyle="round,pad=0.06,rounding_size=0.15",
                       facecolor=color, edgecolor=color, linewidth=1.5,
                       alpha=0.15)
    ax.add_patch(r)
    ax.text(x + w_here/2, y + box_h - 0.25, title,
            ha="center", va="top", fontsize=10, fontweight="bold",
            color=color, family="Helvetica")
    for i, g in enumerate(genes):
        row = i
        ax.text(x + w_here/2, y + box_h - 0.7 - row*0.35, g,
                ha="center", va="center", fontsize=13, fontweight="bold",
                color=NAVY, family="Helvetica")

# arrows / annotation for the outcome
ax.annotate("", xy=(5, 1.5), xytext=(5, 2.2),
            arrowprops=dict(arrowstyle="->", color=NAVY, lw=1.2))
ax.text(5, 1.1,
        "→  hair-follicle neogenesis · scar-free healing",
        ha="center", va="top", fontsize=11, fontweight="bold",
        color=CORAL, family="Helvetica", style="italic")

# citation footer
ax.text(5, 0.15,
        "Signature identified in the Ferreira-lab Mus + Acomys homeostatic "
        "atlas (manuscript in preparation).",
        ha="center", va="bottom", fontsize=8, style="italic", color=MUTED,
        family="Helvetica")

plt.tight_layout()
out = "/Users/apple/Downloads/preview_output/PREVIEW_hero_genes_schema.png"
plt.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
print(f"Wrote {out}")
