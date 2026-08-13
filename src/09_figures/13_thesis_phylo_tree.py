"""
TimeTree-style time-calibrated cladogram — matches the layout the user showed.
9 species, right-angle branches, blue node-age labels, phenotype colour chips.
"""
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib as mpl

mpl.rcParams.update({
    "font.family": "Helvetica",
    "font.size": 10,
    "axes.linewidth": 0.6,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})

INK      = "#111111"
GREY     = "#666666"
NODE_BLU = "#2E60A8"     # TimeTree-style blue for node ages
REGEN    = "#3E8E5F"
SCAR     = "#B84A3E"

# ---- species (tip_y position, top→bottom) --------------------------------
# TimeTree v5 (Kumar et al. 2022)
species = [   # (latin, common, y, phenotype)
    ("Rattus norvegicus",     "Norway rat",         9, SCAR),
    ("Mus musculus",          "house mouse",        8, SCAR),
    ("Acomys cahirinus",      "Cairo spiny mouse",  7, REGEN),
    ("Oryctolagus cuniculus", "European rabbit",    6, REGEN),
    ("Homo sapiens",          "human",              5, SCAR),
    ("Bos taurus",            "cow",                4, SCAR),
    ("Rangifer tarandus",     "reindeer",           3, REGEN),
    ("Ambystoma mexicanum",   "axolotl",            2, REGEN),
    ("Danio rerio",           "zebrafish",          1, REGEN),
]

# ---- internal nodes ------------------------------------------------------
# (age_Mya, y_upper_child, y_lower_child, label_offset_x, label_offset_y)
nodes = [
    dict(age=13,  y_up=9, y_dn=8),  # Rattus + Mus
    dict(age=25,  y_up=8.5, y_dn=7),  # (Rat+Mus) + Acomys
    dict(age=82,  y_up=7.75, y_dn=6),  # + Rabbit
    dict(age=88,  y_up=6.875, y_dn=5),  # + Human
    dict(age=28,  y_up=4, y_dn=3),  # Bos + Reindeer
    dict(age=94,  y_up=5.94, y_dn=3.5),  # Euarchontoglires + Cetartiodactyla
    dict(age=352, y_up=4.72, y_dn=2),  # + Axolotl
    dict(age=429, y_up=3.36, y_dn=1),  # + Zebrafish (root)
]

fig, ax = plt.subplots(figsize=(11.5, 5.2), dpi=200)
ax.set_xlim(475, -235)     # invert; leave room for tip labels on right
ax.set_ylim(0.4, 10.2)

# ---- branches (right-angle cladogram) ------------------------------------
def hline(x0, x1, y, lw=0.9):
    ax.plot([x0, x1], [y, y], color=INK, lw=lw, solid_capstyle="butt", zorder=3)
def vline(x, y0, y1, lw=0.9):
    ax.plot([x, x], [y0, y1], color=INK, lw=lw, solid_capstyle="butt", zorder=3)

# tip → immediate parent x
tip_parent = {9:13, 8:13, 7:25, 6:82, 5:88, 4:28, 3:28, 2:352, 1:429}
for _, _, y, _ in species:
    hline(0, tip_parent[y], y)

# each internal node draws vertical bar + horizontal from prior deeper node
# 13 Mya: Rattus (9) + Mus (8)
vline(13, 9, 8)
# 25 Mya: (Rat/Mus midpoint 8.5) + Acomys (7)
hline(13, 25, 8.5); vline(25, 8.5, 7)
# 82 Mya: prev midpoint (7.75) + Rabbit (6)
hline(25, 82, 7.75); vline(82, 7.75, 6)
# 88 Mya: prev midpoint (6.875) + Human (5)
hline(82, 88, 6.875); vline(88, 6.875, 5)
# 28 Mya (independent): Bos (4) + Reindeer (3)
vline(28, 4, 3)
# 94 Mya: Euarchontoglires (88 → midpoint 5.94) + Cetartiodactyla (28 → 3.5)
hline(88, 94, 5.94); hline(28, 94, 3.5); vline(94, 5.94, 3.5)
# 352 Mya: mammals (94 → 4.72) + Axolotl (2)
hline(94, 352, 4.72); vline(352, 4.72, 2)
# 429 Mya (root): Tetrapoda (352 → 3.36) + Zebrafish (1)
hline(352, 429, 3.36); vline(429, 3.36, 1)

# ---- node-age labels (blue, above node) ----------------------------------
node_labels = [
    (13,  8.5,  "13"),
    (25,  7.75, "25"),
    (82,  6.875,"82"),
    (88,  5.94, "88"),
    (28,  3.5,  "28"),
    (94,  4.72, "94"),
    (352, 3.36, "352"),
    (429, 2.18, "429"),
]
for x, y, txt in node_labels:
    ax.text(x + 4, y + 0.18, txt,
            fontsize=8, color=NODE_BLU, ha="right", va="bottom",
            fontweight="bold")

# ---- root label (subtle, standard practice) -----------------------------
# The 429-Mya root IS the LCA of every species in the tree
# (Osteichthyes = bony vertebrates). Label it minimally, in italics,
# matching TimeTree conventions.
ROOT_X, ROOT_Y = 429, 2.18
ax.scatter(ROOT_X, ROOT_Y, s=28, facecolor="white",
           edgecolor=NODE_BLU, linewidth=1.2, zorder=6)
ax.text(ROOT_X - 8, ROOT_Y, "Osteichthyes",
        fontsize=8, color=NODE_BLU, fontweight="bold",
        style="italic", ha="right", va="center")
ax.text(ROOT_X - 8, ROOT_Y - 0.35, "(last common ancestor)",
        fontsize=7, color=GREY, ha="right", va="center", style="italic")

# ---- tip labels + phenotype chips ---------------------------------------
LABEL_X = -8            # start of Latin name
COMMON_X = -115         # start of common name (fixed column, past longest Latin)
CHIP_X  = -215          # phenotype chip far right (with visible gap)

# column header above the chip column, labelled so reader knows what it is
ax.text(CHIP_X + 3, 10.05, "Skin\nphenotype",
        fontsize=8, fontweight="bold", color=INK,
        ha="center", va="bottom")
for latin, common, y, ph in species:
    ax.text(LABEL_X, y, latin,
            fontsize=10, va="center", ha="left",
            style="italic", color=INK)
    ax.text(COMMON_X, y, f"({common})",
            fontsize=9, va="center", ha="left",
            color=GREY)
    # phenotype chip
    ax.add_patch(Rectangle((CHIP_X, y - 0.22), 6, 0.44,
                           facecolor=ph, edgecolor="black", linewidth=0.5,
                           zorder=5, clip_on=False))

# ---- time axis ----------------------------------------------------------
ax.set_xlabel("Time (million years ago)", fontsize=9.5, color=INK)
ax.set_xticks([0, 50, 100, 150, 200, 250, 300, 350, 400, 450])
ax.tick_params(axis="x", labelsize=8.5, color=GREY, length=3, pad=2)
ax.set_yticks([])
for s in ("top", "right", "left"):
    ax.spines[s].set_visible(False)
ax.spines["bottom"].set_color(INK)
ax.spines["bottom"].set_linewidth(0.6)

# ---- title + subtitle ---------------------------------------------------
fig.text(0.02, 0.96,
    "Phylogenetic relationships among species for the cross-species skin scRNA-seq atlas",
    fontsize=12, fontweight="bold", ha="left", va="top", color=INK)
fig.text(0.02, 0.92,
    "Time-calibrated tree (9 species, all internal nodes labelled)  ·  Source: TimeTree (timetree.org)",
    fontsize=9, ha="left", va="top", color=GREY)

# ---- source citation bottom left ---------------------------------------
fig.text(0.02, 0.02,
    "Kumar S et al. TimeTree 5. Mol Biol Evol 2022;39(8):msac174.",
    fontsize=8, ha="left", va="bottom", color=GREY)

# ---- legend below the plot (horizontal, near x-axis label) --------------
leg_y = 0.045
fig.text(0.55, leg_y + 0.005,
         "Skin phenotype:", fontsize=9, fontweight="bold",
         ha="right", va="center", color=INK)
# regenerator chip
fig.patches.append(Rectangle((0.56, leg_y - 0.007), 0.018, 0.024,
                             transform=fig.transFigure,
                             facecolor=REGEN, edgecolor="black",
                             linewidth=0.5, zorder=11, clip_on=False))
fig.text(0.585, leg_y + 0.005, "Regenerator (scar-free wound healing)",
         fontsize=8.5, va="center", ha="left",
         transform=fig.transFigure, color=INK)
# scarrer chip
fig.patches.append(Rectangle((0.78, leg_y - 0.007), 0.018, 0.024,
                             transform=fig.transFigure,
                             facecolor=SCAR, edgecolor="black",
                             linewidth=0.5, zorder=11, clip_on=False))
fig.text(0.805, leg_y + 0.005, "Scarrer (fibrotic repair)",
         fontsize=8.5, va="center", ha="left",
         transform=fig.transFigure, color=INK)

plt.subplots_adjust(left=0.02, right=0.98, top=0.85, bottom=0.14)
for ext in ("png", "pdf"):
    out = f"/Users/apple/Downloads/preview_output/PREVIEW_phylo_tree_v4.{ext}"
    plt.savefig(out, dpi=600 if ext == "png" else None,
                bbox_inches="tight", facecolor="white")
    print(f"Wrote {out}")
plt.close()
