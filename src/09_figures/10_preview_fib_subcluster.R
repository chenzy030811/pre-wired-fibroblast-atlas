# Fibroblast subclustering — reproduces Sinha 2022 Fig 2C style.
# Loads annotated RDS, subsets fibroblast lineage, re-runs PCA + Harmony + UMAP.
suppressPackageStartupMessages({
  library(Seurat); library(tidyverse); library(patchwork); library(ggrepel)
  library(harmony)
})
options(future.globals.maxSize = 8 * 1024^3)

base_dir <- "/Users/apple/Downloads/研究生毕业论文/reindeer_data"
fib_rds  <- file.path(base_dir, "Reindeer_Day0_Fibroblast_subcluster.rds")

# Fast path -- load the pre-computed subclustered object (saved from prior run)
if (file.exists(fib_rds)) {
  cat("Loading pre-computed fibroblast subset...\n")
  fib <- readRDS(fib_rds)
  # Re-run UMAP -- aim for paper Fig 2C: one connected tree-shaped manifold.
  # Fewer dims + very large n.neighbors + local_connectivity=2 = forces all
  # subclusters to connect to neighbors rather than splitting into islands.
  fib <- RunUMAP(fib, reduction = "harmony", dims = 1:10,
                 reduction.name = "umap",
                 min.dist = 0.3, n.neighbors = 200, spread = 1.0,
                 local.connectivity = 2,
                 verbose = FALSE)
  saveRDS(fib, fib_rds)
} else {
  seu <- readRDS(file.path(base_dir, "Reindeer_Day0_annotated.rds"))
  fib_paper_levels <- intersect(c("Fibroblast", "Fibroblasts"),
                                levels(seu$celltype_paper))
  fib <- subset(seu, celltype_paper %in% fib_paper_levels)
  fib@reductions <- list()
  fib <- JoinLayers(fib)
  fib[["RNA"]] <- split(fib[["RNA"]], f = fib$orig.ident)
  fib <- fib %>%
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(nfeatures = 2000, verbose = FALSE) %>%
    ScaleData(vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE) %>%
    RunPCA(npcs = 30, verbose = FALSE)
  fib <- IntegrateLayers(fib, method = HarmonyIntegration,
                         orig.reduction = "pca", new.reduction = "harmony",
                         verbose = FALSE)
  fib <- fib %>%
    FindNeighbors(reduction = "harmony", dims = 1:15, verbose = FALSE) %>%
    FindClusters(resolution = 0.5, verbose = FALSE) %>%
    RunUMAP(reduction = "harmony", dims = 1:15,
            reduction.name = "umap", verbose = FALSE)
  fib <- JoinLayers(fib)
  saveRDS(fib, fib_rds)
}
cat(sprintf("Fibroblast cells: %d   subclusters: %d\n",
            ncol(fib), length(unique(fib$seurat_clusters))))

# ---- Plot — cluster + parent fine label ----
theme_umap_min <- function(base = 8) {
  theme_classic(base_size = base, base_family = "Helvetica") + theme(
    plot.title  = element_text(size = base + 2, face = "bold", hjust = 0.5),
    axis.text   = element_blank(), axis.ticks = element_blank(),
    axis.title  = element_blank(), axis.line  = element_blank(),
    strip.text  = element_text(size = base + 1, face = "bold"),
    strip.background = element_blank(),
    plot.margin = margin(6, 8, 6, 8)
  )
}
corner_axes <- function(df, frac = 0.10, label_size = 2.4, gap = 0.06) {
  xr <- range(df$UMAP1); yr <- range(df$UMAP2)
  xs <- xr[1] - diff(xr) * gap
  ys <- yr[1] - diff(yr) * gap
  xl <- xs + diff(xr) * frac
  yl <- ys + diff(yr) * frac
  list(
    annotate("segment", x = xs, xend = xl, y = ys, yend = ys,
             arrow = arrow(length = unit(0.08, "cm"), type = "closed"),
             linewidth = 0.4),
    annotate("segment", x = xs, xend = xs, y = ys, yend = yl,
             arrow = arrow(length = unit(0.08, "cm"), type = "closed"),
             linewidth = 0.4),
    annotate("text", x = xs + diff(xr) * frac * 0.5,
             y = ys - diff(yr) * 0.04, label = "UMAP 1",
             size = label_size, hjust = 0.5, vjust = 1),
    annotate("text", x = xs - diff(xr) * 0.03,
             y = ys + diff(yr) * frac * 0.5, label = "UMAP 2",
             size = label_size, hjust = 0.5, vjust = 0.5, angle = 90)
  )
}
expand_low <- function(rng, gap = 0.22) {
  d <- diff(rng); c(rng[1] - d * gap, rng[2])
}

emb <- Embeddings(fib, "umap")

# Panel A: subclusters (numeric labels, paper Fig 2C style)
df_cl <- data.frame(UMAP1 = emb[, 1], UMAP2 = emb[, 2],
                    cluster = fib$seurat_clusters)
centers_cl <- df_cl %>% group_by(cluster) %>%
  summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2), .groups = "drop")

# 12 deep saturated colors (one per fibroblast subcluster)
CLUSTER_COLORS <- c(
  "0"  = "#C66055",  # deep salmon
  "1"  = "#AD5818",  # burnt orange
  "2"  = "#A07C20",  # deep gold
  "3"  = "#7C8E2C",  # deep olive
  "4"  = "#2E7D32",  # deep green
  "5"  = "#2D9686",  # deep teal
  "6"  = "#5090C6",  # sky blue
  "7"  = "#1F4F94",  # deep blue
  "8"  = "#9560AE",  # deep purple
  "9"  = "#B83668",  # rose magenta
  "10" = "#85316C",  # dark magenta
  "11" = "#4A4E69"   # slate
)

# Paper-Fig-2C style metadata: italic "N fibroblasts" tucked above the corner arrow
fib_meta <- sprintf("%s fibroblasts", format(ncol(fib), big.mark = ","))
xr <- range(df_cl$UMAP1); yr <- range(df_cl$UMAP2)
meta_x <- xr[1] - diff(xr) * 0.06
meta_y <- yr[1] + diff(yr) * 0.04   # just above the bottom-left arrow

p_cl <- ggplot(df_cl, aes(UMAP1, UMAP2, color = cluster)) +
  geom_point(size = 1.1, alpha = 0.95, stroke = 0) +
  scale_color_manual(values = CLUSTER_COLORS) +
  ggrepel::geom_text_repel(data = centers_cl, aes(label = cluster),
                           color = "black", size = 5, fontface = "bold",
                           bg.color = "white", bg.r = 0.2,
                           segment.color = NA, max.overlaps = Inf) +
  corner_axes(df_cl) +
  annotate("text", x = meta_x, y = meta_y,
           label = fib_meta, fontface = "italic",
           size = 3.2, hjust = 0, vjust = 0, color = "grey20") +
  coord_cartesian(xlim = expand_low(range(df_cl$UMAP1)),
                  ylim = expand_low(range(df_cl$UMAP2)),
                  clip = "off") +
  theme_umap_min(9) + NoLegend()

ggsave("/tmp/PREVIEW_Fib_subcluster.png",
       plot = p_cl, width = 7, height = 6,
       dpi = 300, units = "in", bg = "white", device = "png", type = "cairo")

# Panel B: colored by parent fine label (shows MP / Fib / Chondro / Tendon)
FIB_FINE_COLORS <- c(
  "Mesenchymal_progenitor"  = "#9560AE",  # purple
  "Fibroblast"              = "#C66055",  # salmon
  "Tendon_fibroblast"       = "#2D9686",  # teal
  "Chondrocyte"             = "#A07C20"   # dark gold
)
df_fine <- data.frame(UMAP1 = emb[, 1], UMAP2 = emb[, 2],
                      fine = fib$celltype_fine)
centers_fine <- df_fine %>% group_by(fine) %>%
  summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2), .groups = "drop")

p_fine <- ggplot(df_fine, aes(UMAP1, UMAP2, color = fine)) +
  geom_point(size = 0.55, alpha = 0.85, stroke = 0) +
  scale_color_manual(values = FIB_FINE_COLORS, name = NULL, drop = TRUE) +
  ggrepel::geom_text_repel(data = centers_fine, aes(label = fine),
                           color = "grey15", size = 2.8,
                           bg.color = "white", bg.r = 0.18,
                           fontface = "plain", segment.color = NA,
                           max.overlaps = Inf) +
  corner_axes(df_fine) +
  coord_cartesian(xlim = expand_low(range(df_fine$UMAP1)),
                  ylim = expand_low(range(df_fine$UMAP2)),
                  clip = "off") +
  labs(title = "Fibroblast subtypes (parent fine labels)") +
  theme_umap_min(9) + NoLegend()

ggsave("/tmp/PREVIEW_Fib_subcluster_fine.png",
       plot = p_fine, width = 7, height = 6,
       dpi = 300, units = "in", bg = "white", device = "png", type = "cairo")

cat("\nDone. Wrote:\n",
    "/tmp/PREVIEW_Fib_subcluster.png\n",
    "/tmp/PREVIEW_Fib_subcluster_fine.png\n")
