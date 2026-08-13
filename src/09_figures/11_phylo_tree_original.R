# ============================================================
# Phylogenetic tree for cross-species skin scRNA-seq atlas
# 9 species, 8 internal nodes labeled, time-calibrated.
# Divergence times from TimeTree.org (Kumar et al. 2022).
# ============================================================
suppressPackageStartupMessages({
  library(ape); library(ggtree); library(ggplot2); library(dplyr); library(tibble)
})

# ---- 1. Time-calibrated Newick tree (branch lengths in Mya) ----
newick <- "(Danio_rerio:429,(Ambystoma_mexicanum:352,((Rangifer_tarandus:28,Bos_taurus:28):66,((((Mus_musculus:13,Rattus_norvegicus:13):12,Acomys_cahirinus:25):57,Oryctolagus_cuniculus:82):6,Homo_sapiens:88):6):258):77);"
tree <- read.tree(text = newick)

# ---- 2. Tip annotation ----
tip_info <- tibble(
  label   = c("Danio_rerio", "Ambystoma_mexicanum",
              "Rangifer_tarandus", "Bos_taurus", "Homo_sapiens",
              "Oryctolagus_cuniculus", "Acomys_cahirinus",
              "Mus_musculus", "Rattus_norvegicus"),
  common  = c("zebrafish", "axolotl", "reindeer", "cow", "human",
              "European rabbit", "Cairo spiny mouse",
              "house mouse", "Norway rat"),
  pheno   = c("Regenerator", "Regenerator", "Regenerator", "Scarrer", "Scarrer",
              "Regenerator", "Regenerator", "Scarrer", "Scarrer")
)

# ---- 3. Internal node ages ----
n_tip   <- length(tree$tip.label)
n_node  <- tree$Nnode
node_id <- (n_tip + 1):(n_tip + n_node)
heights <- node.depth.edgelength(tree)
max_h   <- max(heights[1:n_tip])
node_age <- max_h - heights[node_id]
tree$node.label <- as.character(round(node_age, 0))

# ---- 4. Palette ----
PHENO_COLORS <- c(Regenerator = "#2E7D32", Scarrer = "#C0392B")

# ---- 5. Plot (ggtree) ----
# Attach tip metadata to the tree; revts() flips so present = 0 on right
p <- revts(ggtree(tree, size = 0.55, color = "grey15") %<+% tip_info)

# Tip labels: Latin (italic) + (common name)  -- simple, no aes parsing
p <- p +
  geom_tiplab(aes(label = paste0(gsub("_", " ", label), "  (", common, ")")),
              hjust = 0, offset = 4, size = 3.3,
              color = "grey15", fontface = "italic") +
  # Phenotype square to the right of each tip (after revts(), x = 0 at present)
  geom_tippoint(aes(color = pheno), shape = 15, size = 4,
                position = position_nudge(x = 145)) +
  # Tile column extends past 0; coord_cartesian below will keep it visible
  scale_color_manual(values = PHENO_COLORS, name = "Skin phenotype",
                     guide = guide_legend(override.aes = list(size = 4))) +
  # Internal-node divergence-time labels (stored in tree$node.label)
  geom_nodelab(size = 3.2, hjust = 1.15, vjust = -0.5,
               color = "#1F4F94", fontface = "bold") +
  # X-axis = Mya (after revts() values are negative; relabel as positive)
  scale_x_continuous(
    breaks = -seq(0, 450, by = 50),
    labels = seq(0, 450, by = 50),
    name   = "Time (million years ago)") +
  coord_cartesian(xlim = c(-470, 165), clip = "off") +
  labs(title    = "Phylogenetic relationships among species for the cross-species skin scRNA-seq atlas",
       subtitle = "Time-calibrated tree (9 species, all internal nodes labelled)  -  Source: TimeTree (timetree.org)",
       caption  = "Kumar S et al. TimeTree 5. Mol Biol Evol 2022;39(8):msac174.") +
  theme_tree2(base_size = 9, base_family = "Helvetica") +
  theme(
    plot.title       = element_text(size = 12, face = "bold", hjust = 0),
    plot.subtitle    = element_text(size = 9,  color = "grey30", hjust = 0,
                                    margin = margin(b = 10)),
    plot.caption     = element_text(size = 8,  color = "grey50", hjust = 0),
    axis.title.x     = element_text(size = 10, margin = margin(t = 4)),
    axis.text.x      = element_text(size = 9,  color = "black"),
    axis.line.x      = element_line(linewidth = 0.4, color = "black"),
    axis.ticks.x     = element_line(linewidth = 0.3, color = "black"),
    legend.position  = c(0.92, 0.20),
    legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.2),
    legend.title     = element_text(size = 9, face = "bold"),
    legend.text      = element_text(size = 8),
    plot.margin      = margin(10, 14, 8, 12)
  )

ggsave("/tmp/PREVIEW_phylo_tree.png", plot = p,
       width = 13, height = 5.5, dpi = 300, units = "in",
       bg = "white", device = "png", type = "cairo")
cat("Wrote /tmp/PREVIEW_phylo_tree.png\n")
