# Quick preview of pastel UMAP styling — loads annotated RDS, no recompute.
suppressPackageStartupMessages({
  library(Seurat); library(tidyverse); library(patchwork); library(ggrepel)
})

base_dir <- "/Users/apple/Downloads/研究生毕业论文/reindeer_data"
seu <- readRDS(file.path(base_dir, "Reindeer_Day0_annotated.rds"))

CELLTYPE_COLORS <- c(
  # Deeper than Fig 2A — punchier for thesis/PPT
  "Fibroblast"              = "#C66055",
  "Fibroblasts"             = "#C66055",
  "Basal Keratinocyte"      = "#A07C20",
  "Basal Keratinocytes"     = "#A07C20",
  "Suprabasal Keratinocyte" = "#AD5818",
  "Endothelium"             = "#2D9686",
  "Endothelial"             = "#2D9686",
  "Melanocyte"              = "#7C8E2C",
  "T cell"                  = "#9560AE",
  "T-cell"                  = "#9560AE",
  "CD45+IL1β+ Myeloid"      = "#5090C6",
  "Myeloid"                 = "#5090C6",
  "VSM"                     = "#1F4F94",
  "Macrophage (MΦ)"         = "#B83668",
  "Macrophages"             = "#B83668",
  "Schwann"                 = "#85316C"
)
TISSUE_COLORS <- c(Velvet = "#E07A5F", Back = "#3D8DAE")

theme_umap <- function(base = 8) {
  theme_classic(base_size = base, base_family = "Helvetica") + theme(
    plot.title  = element_text(size = base + 2, face = "bold", hjust = 0.5),
    axis.text   = element_blank(),
    axis.ticks  = element_blank(),
    axis.title  = element_blank(),
    axis.line   = element_blank(),
    strip.text  = element_text(size = base + 1, face = "bold"),
    strip.background = element_blank(),
    plot.margin = margin(6, 8, 6, 8)
  )
}

# Bottom-left UMAP 1 / UMAP 2 corner arrow — placed BELOW/LEFT of all data
# Caller must use coord_cartesian with expand_low() so the arrows fit.
corner_axes <- function(df, frac = 0.10, label_size = 2.4, gap = 0.06) {
  xr <- range(df$UMAP1); yr <- range(df$UMAP2)
  xs <- xr[1] - diff(xr) * gap
  ys <- yr[1] - diff(yr) * gap
  xl <- xs + diff(xr) * frac
  yl <- ys + diff(yr) * frac
  list(
    annotate("segment", x = xs, xend = xl, y = ys, yend = ys,
             arrow = arrow(length = unit(0.08, "cm"), type = "closed"),
             color = "black", linewidth = 0.4),
    annotate("segment", x = xs, xend = xs, y = ys, yend = yl,
             arrow = arrow(length = unit(0.08, "cm"), type = "closed"),
             color = "black", linewidth = 0.4),
    annotate("text", x = xs + diff(xr) * frac * 0.5,
             y = ys - diff(yr) * 0.04,
             label = "UMAP 1", size = label_size, hjust = 0.5, vjust = 1),
    annotate("text", x = xs - diff(xr) * 0.03,
             y = ys + diff(yr) * frac * 0.5,
             label = "UMAP 2", size = label_size, hjust = 0.5, vjust = 0.5,
             angle = 90)
  )
}

# Helper — expand coord range on low side to make room for corner_axes + labels
expand_low <- function(rng, gap = 0.22) {
  d <- diff(rng); c(rng[1] - d * gap, rng[2])
}

# Ensure paper-aligned label column exists
if (!"celltype_paper" %in% colnames(seu@meta.data)) {
  stop("celltype_paper not found -- please run the full pipeline once first.")
}
Idents(seu) <- seu$celltype_paper

paper_palette <- CELLTYPE_COLORS[
  intersect(levels(seu$celltype_paper), names(CELLTYPE_COLORS))]
miss <- setdiff(levels(seu$celltype_paper), names(paper_palette))
if (length(miss)) paper_palette <- c(paper_palette,
                                     setNames(rep("#888888", length(miss)), miss))

# ---- Single-panel pastel UMAP ----
emb <- Embeddings(seu, "umap")
df  <- data.frame(UMAP1 = emb[, 1], UMAP2 = emb[, 2],
                  grp   = seu$celltype_paper)
centers <- df %>% group_by(grp) %>%
  summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2), .groups = "drop")

p_single <- ggplot(df, aes(UMAP1, UMAP2, color = grp)) +
  geom_point(size = 0.6, alpha = 0.85, stroke = 0) +
  scale_color_manual(values = paper_palette, name = NULL, drop = FALSE) +
  ggrepel::geom_text_repel(data = centers, aes(label = grp),
                           color = "grey15", size = 2.8,
                           bg.color = "white", bg.r = 0.18,
                           fontface = "plain", segment.color = NA,
                           max.overlaps = Inf) +
  corner_axes(df) +
  coord_cartesian(xlim = expand_low(range(df$UMAP1)),
                  ylim = expand_low(range(df$UMAP2)),
                  clip = "off") +
  labs(title = "Cell types - paper-aligned (preview)") +
  theme_umap(8) + NoLegend()

ggsave("/tmp/PREVIEW_UMAP_paper.png",
       plot = p_single, width = 7.5, height = 6,
       dpi = 300, units = "in", bg = "white", device = "png", type = "cairo")

# ---- Same plot, NO text labels (for hand-annotation) ----
p_single_blank <- ggplot(df, aes(UMAP1, UMAP2, color = grp)) +
  geom_point(size = 0.6, alpha = 0.85, stroke = 0) +
  scale_color_manual(values = paper_palette, name = NULL, drop = FALSE) +
  corner_axes(df) +
  coord_cartesian(xlim = expand_low(range(df$UMAP1)),
                  ylim = expand_low(range(df$UMAP2)),
                  clip = "off") +
  labs(title = "Cell types - paper-aligned (unlabeled)") +
  theme_umap(8) + NoLegend()
ggsave("/tmp/PREVIEW_UMAP_paper_blank.png",
       plot = p_single_blank, width = 7.5, height = 6,
       dpi = 300, units = "in", bg = "white", device = "png", type = "cairo")

# ---- Split panel (Back | Velvet) ----
seu$tissue_paper <- factor(
  ifelse(seu$orig.ident == "Day0_Antler", "Velvet", "Back"),
  levels = c("Back", "Velvet"))

xlims <- expand_low(range(emb[, 1]) + c(-0.5, 0.5))
ylims <- expand_low(range(emb[, 2]) + c(-0.5, 0.5))

mk <- function(tissue, border) {
  d <- df; d$tissue_paper <- seu$tissue_paper
  d <- d[d$tissue_paper == tissue, ]
  cc <- d %>% group_by(grp) %>%
    summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2), .groups = "drop")
  ggplot(d, aes(UMAP1, UMAP2, color = grp)) +
    geom_point(size = 0.6, alpha = 0.85, stroke = 0) +
    scale_color_manual(values = paper_palette, drop = FALSE) +
    ggrepel::geom_text_repel(data = cc, aes(label = grp),
                             color = "grey15", size = 2.9,
                             bg.color = "white", bg.r = 0.18,
                             fontface = "plain", segment.color = NA,
                             max.overlaps = Inf) +
    corner_axes(data.frame(UMAP1 = range(emb[, 1]),
                            UMAP2 = range(emb[, 2]))) +
    coord_cartesian(xlim = xlims, ylim = ylims, clip = "off") +
    labs(title = tissue) +
    NoLegend() + theme_umap(9) +
    theme(plot.title   = element_text(face = "bold", size = 13,
                                      hjust = 0.02, color = border),
          panel.border = element_rect(color = border, fill = NA,
                                      linewidth = 1.2))
}

p_split <- (mk("Back", TISSUE_COLORS[["Back"]]) |
            mk("Velvet", TISSUE_COLORS[["Velvet"]])) +
  plot_annotation(
    title = "Day 0 - paper-aligned (preview, pastel)",
    theme = theme(plot.title = element_text(size = 12, face = "bold", hjust = 0))) &
  theme(legend.position = "none")

ggsave("/tmp/PREVIEW_UMAP_paper_split.png",
       plot = p_split, width = 11, height = 5.5,
       dpi = 300, units = "in", bg = "white", device = "png", type = "cairo")

# ---- Split, NO text labels ----
mk_blank <- function(tissue, border) {
  d <- df; d$tissue_paper <- seu$tissue_paper
  d <- d[d$tissue_paper == tissue, ]
  ggplot(d, aes(UMAP1, UMAP2, color = grp)) +
    geom_point(size = 0.6, alpha = 0.85, stroke = 0) +
    scale_color_manual(values = paper_palette, drop = FALSE) +
    corner_axes(data.frame(UMAP1 = range(emb[, 1]),
                            UMAP2 = range(emb[, 2]))) +
    coord_cartesian(xlim = xlims, ylim = ylims, clip = "off") +
    labs(title = tissue) +
    NoLegend() + theme_umap(9) +
    theme(plot.title   = element_text(face = "bold", size = 13,
                                      hjust = 0.02, color = border),
          panel.border = element_rect(color = border, fill = NA,
                                      linewidth = 1.2))
}
p_split_blank <- (mk_blank("Back", TISSUE_COLORS[["Back"]]) |
                  mk_blank("Velvet", TISSUE_COLORS[["Velvet"]])) +
  plot_annotation(
    title = "Day 0 - paper-aligned (unlabeled, pastel)",
    theme = theme(plot.title = element_text(size = 12, face = "bold", hjust = 0))) &
  theme(legend.position = "none")
ggsave("/tmp/PREVIEW_UMAP_paper_split_blank.png",
       plot = p_split_blank, width = 11, height = 5.5,
       dpi = 300, units = "in", bg = "white", device = "png", type = "cairo")

cat("Done. Wrote:\n",
    "/tmp/PREVIEW_UMAP_paper.png", "\n",
    "/tmp/PREVIEW_UMAP_paper_split.png", "\n")
