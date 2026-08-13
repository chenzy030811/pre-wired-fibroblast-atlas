# Reindeer Day 0 scRNA-seq pipeline 
# Data    : GSE168748, Sinha et al. 2022 Cell (14,540 cells total)
# Samples : Day0_Antler (velvet, regenerative) + Day0_Back (dorsal, scarring)
# Author  : Zhiyi Chen

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(harmony)        # needed for HarmonyIntegration
})

# Raise future's per-worker memory cap (default 500 MiB is too small for
# Seurat v5 IntegrateLayers on ~10k cells with full-gene scale.data).
options(future.globals.maxSize = 8 * 1024^3)   # 8 GB


# 0. Global config (single source of truth)
base_dir <- "/Users/apple/Downloads/研究生毕业论文/reindeer_data"
samples  <- c("Day0_Antler", "Day0_Back")

# QC thresholds. Lower bounds prevent low-feature cells from forming
# an artefact cluster downstream.(基因数、UMI 总数和线粒体基因比例)
TH_FEATURE_LOW  <- 800
TH_FEATURE_HIGH <- 4500
TH_COUNT_LOW    <- 1500
TH_COUNT_HIGH   <- 15000
TH_MT           <- 5

# Reindeer is annotated against the bovine reference, so mitochondrial
# genes use gene-symbol naming (ND/COX/ATP/CYTB) rather than "MT-" prefix.
MT_PATTERN <- "^MT-|^mt-|^Mt-|^ND[1-6]$|^ND4L$|^COX[1-3]$|^ATP[68]$|^CYTB$"

SAMPLE_COLORS <- c(Day0_Antler = "#E07A5F",
                   Day0_Back   = "#3D8DAE")
TISSUE_COLORS <- c(Velvet = "#E07A5F", Back = "#3D8DAE")


# ============================================================
# Publication design system (Cell-journal aesthetic)
# Color, theme and save helpers consumed by every figure below.
# ============================================================
suppressPackageStartupMessages({
  library(ggrepel)
  library(scales)
})

# Cell-type palette — deep saturated (final approved)
CELLTYPE_COLORS <- c(
  "Fibroblast"              = "#C66055",
  "Basal Keratinocyte"      = "#A07C20",
  "Suprabasal Keratinocyte" = "#AD5818",
  "Endothelium"             = "#2D9686",
  "Melanocyte"              = "#7C8E2C",
  "T cell"                  = "#9560AE",
  "CD45+IL1β+ Myeloid" = "#5090C6",
  "VSM"                     = "#1F4F94",
  "Macrophage (MΦ)"    = "#B83668",
  "Schwann"                 = "#85316C",
  # Saved-RDS plural variants -> same hue
  "Fibroblasts"             = "#C66055",
  "Basal Keratinocytes"     = "#A07C20",
  "Endothelial"             = "#2D9686",
  "T-cell"                  = "#9560AE",
  "Myeloid"                 = "#5090C6",
  "Macrophages"             = "#B83668"
)

# Fine-label palette (sub-types share parent hue, deep saturated)
FINE_COLORS <- c(
  "Fibroblast"              = "#C66055",
  "Mesenchymal_progenitor"  = "#A04035",
  "Tendon_fibroblast"       = "#E07A4E",
  "Chondrocyte"             = "#8B3A2F",
  "Basal_keratinocyte"      = "#A07C20",
  "Hair_shaft_keratinocyte" = "#7E601A",
  "Suprabasal_keratinocyte" = "#AD5818",
  "Sebaceous_cell"          = "#7A4818",
  "Endothelial"             = "#2D9686",
  "Lymphatic_EC"            = "#1F7868",
  "Pericyte_SMC"            = "#1F4F94",
  "T_cell"                  = "#9560AE",
  "Dendritic_cell"          = "#5090C6",
  "Macrophage"              = "#B83668",
  "Mast_cell"               = "#D44B82",
  "Melanocyte"              = "#7C8E2C",
  "Schwann"                 = "#85316C"
)

# Publication base theme
theme_pub <- function(base = 8) {
  theme_classic(base_size = base, base_family = "Helvetica") +
    theme(
      plot.title       = element_text(size = base + 2, face = "bold",
                                      hjust = 0.5, margin = margin(b = 4)),
      plot.subtitle    = element_text(size = base, color = "grey30",
                                      hjust = 0.5, margin = margin(b = 6)),
      axis.title       = element_text(size = base + 1),
      axis.text        = element_text(size = base, color = "black"),
      axis.line        = element_line(linewidth = 0.4, color = "black"),
      axis.ticks       = element_line(linewidth = 0.3, color = "black"),
      legend.title     = element_text(size = base, face = "bold"),
      legend.text      = element_text(size = base - 1),
      legend.key.size  = unit(0.35, "cm"),
      strip.text       = element_text(size = base + 1, face = "bold"),
      strip.background = element_blank(),
      plot.margin      = margin(6, 8, 6, 8)
    )
}

# UMAP variant — keeps axis titles but hides ticks/text
theme_umap <- function(base = 8) {
  theme_pub(base) + theme(
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(size = base, color = "grey40")
  )
}

# UMAP minimal — paper Fig 2A style: NO axis lines/ticks/titles at all
# (caller is expected to add corner_axes() for the bottom-left UMAP1/UMAP2 arrow)
theme_umap_minimal <- function(base = 8) {
  theme_pub(base) + theme(
    axis.text   = element_blank(),
    axis.ticks  = element_blank(),
    axis.title  = element_blank(),
    axis.line   = element_blank(),
    plot.margin = margin(8, 10, 8, 10)
  )
}

# Bottom-left corner UMAP 1 / UMAP 2 arrow indicator (paper Fig 2A style).
# Returns a list of ggplot layers; caller does `p + corner_axes(p_data)`.
# `df` must have UMAP1 and UMAP2 columns; we infer plot range from them.
corner_axes <- function(df, frac = 0.18, label_size = 2.4) {
  xr <- range(df$UMAP1); yr <- range(df$UMAP2)
  xs <- xr[1] - diff(xr) * 0.04          # x start (left of data)
  ys <- yr[1] - diff(yr) * 0.04          # y start (below data)
  xl <- xs + diff(xr) * frac             # arrow length
  yl <- ys + diff(yr) * frac
  list(
    annotate("segment", x = xs, xend = xl, y = ys, yend = ys,
             arrow = arrow(length = unit(0.08, "cm"), type = "closed"),
             color = "black", linewidth = 0.4),
    annotate("segment", x = xs, xend = xs, y = ys, yend = yl,
             arrow = arrow(length = unit(0.08, "cm"), type = "closed"),
             color = "black", linewidth = 0.4),
    annotate("text", x = xs + diff(xr) * frac * 0.45, y = ys - diff(yr) * 0.025,
             label = "UMAP 1", size = label_size, hjust = 0.5, vjust = 1,
             color = "black"),
    annotate("text", x = xs - diff(xr) * 0.025, y = ys + diff(yr) * frac * 0.45,
             label = "UMAP 2", size = label_size, hjust = 0.5, vjust = 0.5,
             color = "black", angle = 90)
  )
}

# Save helper -- PNG @ 600 dpi, white background (Word-friendly)
save_pub <- function(plot, name, w, h) {
  ggsave(file.path(base_dir, paste0(name, ".png")),
         plot = plot, width = w, height = h,
         dpi = 600, units = "in", bg = "white")
  invisible(NULL)
}


# 1. Load samples, light pre-filter, merge
seu_list <- lapply(samples, function(s) {
  counts <- Read10X(data.dir = file.path(base_dir, s))
  CreateSeuratObject(counts, project = s,
                     min.cells = 3, min.features = 200)
})

seu_pre <- merge(seu_list[[1]], seu_list[[2]],
                 add.cell.ids = samples,
                 project = "Reindeer_Day0")

seu_pre$tissue <- ifelse(grepl("Antler", seu_pre$orig.ident),
                         "Antler", "Back")
seu_pre[["percent.mt"]] <- PercentageFeatureSet(seu_pre, pattern = MT_PATTERN)


# 2. Apply QC filter
seu_post <- subset(seu_pre,
                   subset = nFeature_RNA > TH_FEATURE_LOW &
                            nFeature_RNA < TH_FEATURE_HIGH &
                            nCount_RNA   > TH_COUNT_LOW   &
                            nCount_RNA   < TH_COUNT_HIGH  &
                            percent.mt   < TH_MT)

cat("Antler:", sum(seu_post$orig.ident == "Day0_Antler"), "\n")
cat("Back:  ", sum(seu_post$orig.ident == "Day0_Back"),   "\n")
cat("Total: ", ncol(seu_post), "/", ncol(seu_pre),
    sprintf("(%.1f%%)\n", 100 * ncol(seu_post) / ncol(seu_pre)))

saveRDS(seu_post, file.path(base_dir, "Reindeer_Day0_postQC.rds"))


# 3. QC figures (publication-grade, Cell aesthetic)

# Publication QC violin: samples side-by-side, dashed thresholds, no jitter,
# Y-axis trimmed so the threshold region dominates the visual.
qc_violin <- function(df, metric, ylab, low = NULL, high = NULL,
                      y_cap = NULL, title = NULL) {
  p <- ggplot(df, aes(x = orig.ident, y = .data[[metric]], fill = orig.ident)) +
    geom_violin(scale = "width", width = 0.85,
                color = "grey25", linewidth = 0.35, alpha = 0.85) +
    geom_boxplot(width = 0.12, outlier.shape = NA,
                 color = "grey15", fill = "white", linewidth = 0.3) +
    scale_fill_manual(values = SAMPLE_COLORS, guide = "none") +
    scale_x_discrete(labels = c(Day0_Antler = "Velvet", Day0_Back = "Back")) +
    labs(x = NULL, y = ylab, title = title) +
    theme_pub(8) +
    theme(axis.title.x = element_blank())
  if (!is.null(low))   p <- p + geom_hline(yintercept = low,
                                           linetype = "dashed",
                                           color = "grey25", linewidth = 0.35)
  if (!is.null(high))  p <- p + geom_hline(yintercept = high,
                                           linetype = "dashed",
                                           color = "grey25", linewidth = 0.35)
  if (!is.null(y_cap)) p <- p + coord_cartesian(ylim = c(0, y_cap))
  p
}

pre_df  <- seu_pre@meta.data
post_df <- seu_post@meta.data

# Y caps so the threshold region is visible (avoid tail-stretched plots)
ycap_feature <- quantile(pre_df$nFeature_RNA, 0.995)
ycap_count   <- quantile(pre_df$nCount_RNA,   0.995)
ycap_mt      <- max(15, TH_MT * 2)

# Before-QC row
p1a <- qc_violin(pre_df, "nFeature_RNA", "Genes / cell",
                 low = TH_FEATURE_LOW, high = TH_FEATURE_HIGH,
                 y_cap = ycap_feature)
p1b <- qc_violin(pre_df, "nCount_RNA",   "UMIs / cell",
                 low = TH_COUNT_LOW,   high = TH_COUNT_HIGH,
                 y_cap = ycap_count)
p1c <- qc_violin(pre_df, "percent.mt",   "Mitochondrial %",
                 high = TH_MT, y_cap = ycap_mt)

# After-QC row (no thresholds; just show that distributions are clean)
p2a <- qc_violin(post_df, "nFeature_RNA", "Genes / cell")
p2b <- qc_violin(post_df, "nCount_RNA",   "UMIs / cell")
p2c <- qc_violin(post_df, "percent.mt",   "Mitochondrial %")

# Compact header (table-style, one row each)
hdr <- sprintf(
  "Day 0 reindeer scRNA-seq | thresholds: nFeature %d-%d, nCount %d-%d, %%mt < %d%%
Before QC: %d cells (Velvet %d / Back %d)   After QC: %d cells (Velvet %d / Back %d, %.1f%% retained)",
  TH_FEATURE_LOW, TH_FEATURE_HIGH, TH_COUNT_LOW, TH_COUNT_HIGH, TH_MT,
  ncol(seu_pre),  sum(seu_pre$orig.ident == "Day0_Antler"),
  sum(seu_pre$orig.ident == "Day0_Back"),
  ncol(seu_post), sum(seu_post$orig.ident == "Day0_Antler"),
  sum(seu_post$orig.ident == "Day0_Back"),
  100 * ncol(seu_post) / ncol(seu_pre)
)

# Tag leftmost plot of each row with row label, baked into the title
p1a <- p1a + labs(title = "Before QC")
p2a <- p2a + labs(title = "After QC")

fig_qc <- (p1a | p1b | p1c) / (p2a | p2b | p2c) +
  plot_annotation(
    title    = "Day 0 reindeer scRNA-seq — quality control",
    subtitle = hdr,
    theme = theme(
      plot.title    = element_text(size = 11, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 7.5, color = "grey30", hjust = 0,
                                   lineheight = 1.3, family = "Helvetica")
    )
  ) +
  plot_layout(heights = c(1, 1))

save_pub(fig_qc, "QC_report_full", w = 10, h = 6.5)


# Per-sample QC violin (kept-cells highlight, single-color, clean).
# Cleaner alternative to the messy 2-color overlay: shows distribution of all
# pre-QC cells with dashed thresholds; subtitle reports kept fraction.
make_qc_violin_clean <- function(df_pre, sample_name, metric,
                                 low = NULL, high = NULL, ylab) {
  d <- df_pre[df_pre$orig.ident == sample_name, ]
  color <- SAMPLE_COLORS[[sample_name]]

  kept <- rep(TRUE, nrow(d))
  if (!is.null(low))  kept <- kept & d[[metric]] >= low
  if (!is.null(high)) kept <- kept & d[[metric]] <= high
  pct_kept <- 100 * sum(kept) / nrow(d)

  y_cap <- if (!is.null(high)) high * 1.15
           else as.numeric(quantile(d[[metric]], 0.995, na.rm = TRUE))

  tissue_label <- ifelse(sample_name == "Day0_Antler", "Velvet", "Back")
  p <- ggplot(d, aes(x = orig.ident, y = .data[[metric]])) +
    geom_violin(fill = color, alpha = 0.85, color = "grey20",
                linewidth = 0.35, scale = "width", width = 0.85) +
    geom_boxplot(width = 0.12, outlier.shape = NA,
                 color = "grey15", fill = "white", linewidth = 0.3) +
    coord_cartesian(ylim = c(0, y_cap)) +
    labs(title    = tissue_label,
         subtitle = sprintf("kept %d / %d  (%.1f%%)",
                            sum(kept), nrow(d), pct_kept),
         x = NULL, y = ylab) +
    theme_pub(8) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  if (!is.null(low))  p <- p + geom_hline(yintercept = low,
                                          linetype = "dashed",
                                          color = "grey25", linewidth = 0.35)
  if (!is.null(high)) p <- p + geom_hline(yintercept = high,
                                          linetype = "dashed",
                                          color = "grey25", linewidth = 0.35)
  p
}

metric_specs <- list(
  list(name = "nFeature_RNA", low = TH_FEATURE_LOW, high = TH_FEATURE_HIGH, ylab = "Genes / cell"),
  list(name = "nCount_RNA",   low = TH_COUNT_LOW,   high = TH_COUNT_HIGH,   ylab = "UMIs / cell"),
  list(name = "percent.mt",   low = NULL,           high = TH_MT,           ylab = "Mitochondrial %")
)

# 2 samples x 3 metrics grid
qc_plots <- list()
for (s in samples) {
  for (m in metric_specs) {
    qc_plots[[length(qc_plots) + 1]] <-
      make_qc_violin_clean(pre_df, s, m$name, m$low, m$high, m$ylab)
  }
}

fig_qc_persample <- wrap_plots(qc_plots, ncol = 3) +
  plot_annotation(
    title = "Per-sample QC distributions (before filtering, with thresholds)",
    theme = theme(plot.title = element_text(size = 11, face = "bold", hjust = 0))
  )

save_pub(fig_qc_persample, "QC_violins_combined", w = 9, h = 7)

# QC result: 9,412 / 14,540 cells (64.7%) retained -> 3,175 Antler + 6,237 Back.
# Violin plots show clean central peaks for nCount and percent.mt; nFeature has
# a minor lower shoulder from a dense granulocyte cluster (data-intrinsic).


# 4. Normalization, HVG selection, scaling
# LogNormalize (scale factor 10,000) corrects for sequencing depth.
# 2,000 HVGs selected by vst, then z-scored across ALL genes with
# nCount_RNA + percent.mt regressed out to remove technical effects.

seu <- seu_post
cat("Loaded object:", ncol(seu), "cells x", nrow(seu), "genes\n")
cat("Samples:\n"); print(table(seu$orig.ident))

seu <- NormalizeData(seu,
                     normalization.method = "LogNormalize",
                     scale.factor         = 1e4)

seu <- FindVariableFeatures(seu,
                            selection.method = "vst",
                            nfeatures        = 2000)

top10 <- head(VariableFeatures(seu), 10)
cat("\nTop 10 HVGs:\n"); print(top10)

# FIGURE 3 - HVG mean-variance plot (publication-grade)
# Build from HVF metadata directly so we can control colors/labels cleanly.
hvf_df <- HVFInfo(seu) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::mutate(
    variable = gene %in% VariableFeatures(seu),
    label    = ifelse(gene %in% top10, gene, NA)
  ) %>%
  dplyr::filter(mean > 0)   # drop genes with avg expr = 0 (log undefined)

n_var    <- sum(hvf_df$variable)
n_nonvar <- sum(!hvf_df$variable)

p_hvg <- ggplot(hvf_df, aes(x = mean, y = variance.standardized,
                            color = variable)) +
  geom_point(size = 0.5, alpha = 0.8) +
  scale_color_manual(values = c(`TRUE` = "#E07A5F", `FALSE` = "grey55"),
                     labels = c(`TRUE`  = sprintf("Variable (n = %d)", n_var),
                                `FALSE` = sprintf("Non-variable (n = %d)", n_nonvar)),
                     breaks = c("TRUE", "FALSE"),
                     name = NULL) +
  scale_x_log10(labels = label_log(),
                breaks = scales::trans_breaks("log10", function(x) 10^x)) +
  ggrepel::geom_text_repel(aes(label = label),
                           color = "black", size = 2.6,
                           bg.color = "white", bg.r = 0.15,
                           min.segment.length = 0, max.overlaps = 20,
                           segment.color = "grey40", segment.size = 0.25,
                           na.rm = TRUE) +
  labs(x = "Average expression (log10)", y = "Standardized variance",
       title = "Highly variable genes (vst, top 2,000)") +
  theme_pub(8) +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = alpha("white", 0.8),
                                         color = NA))

suppressWarnings(save_pub(p_hvg, "HVG_plot", w = 7, h = 4.5))

# Scale ALL genes (so non-HVG markers can be plotted later);
# regress out library size and mitochondrial percentage.
seu <- ScaleData(seu,
                 features        = rownames(seu),
                 vars.to.regress = c("nCount_RNA", "percent.mt"))


# 5. PCA
# Show clean biological signals on PC1-5 (epithelial, fibroblast, immune,
# endothelial, mast cells).

seu <- RunPCA(seu, npcs = 30)

# FIGURE 4 - Elbow plot for choosing # of PCs
sdev   <- Stdev(seu, "pca")[1:30]
elbow_df <- data.frame(PC = seq_along(sdev), SD = sdev)
p_elbow <- ggplot(elbow_df, aes(PC, SD)) +
  geom_line(color = "grey60", linewidth = 0.4) +
  geom_point(color = "#3D8DAE", size = 1.6) +
  geom_vline(xintercept = 10, linetype = "dashed",
             color = "#E07A5F", linewidth = 0.4) +
  annotate("text", x = 10.5, y = max(sdev) * 0.95,
           label = "Selected: 10 PCs", hjust = 0, size = 2.8,
           color = "#E07A5F", fontface = "bold") +
  scale_x_continuous(breaks = seq(0, 30, 5)) +
  labs(x = "Principal component", y = "Standard deviation",
       title = "PCA elbow") +
  theme_pub(8)
save_pub(p_elbow, "ElbowPlot", w = 5, h = 3.5)

# FIGURE 5 - PCA colored by sample (pre-integration batch view)
p_pca <- DimPlot(seu, reduction = "pca", group.by = "orig.ident",
                 raster = TRUE, pt.size = 0.4) +
  scale_color_manual(values = SAMPLE_COLORS,
                     labels = c(Day0_Antler = "Velvet", Day0_Back = "Back"),
                     name = NULL) +
  labs(title = "PCA (no integration)", x = "PC 1", y = "PC 2") +
  theme_pub(8) +
  theme(legend.position = c(0.98, 0.02),
        legend.justification = c(1, 0))
save_pub(p_pca, "PCA_by_sample", w = 5, h = 4.5)


# 5.5 Integration (Harmony / RPCA / CCA comparison)
# Even with 2 samples, integration aligns shared cell types across libraries
# and removes residual technical variation between Antler and Back runs.
# We run all three methods, compare UMAPs, then use Harmony as the default
# (fast, robust, preserves antler-specific populations like ABPCs).

seu <- IntegrateLayers(seu, method = HarmonyIntegration,
                       orig.reduction = "pca",
                       new.reduction  = "harmony",
                       verbose = FALSE)

seu <- IntegrateLayers(seu, method = RPCAIntegration,
                       orig.reduction = "pca",
                       new.reduction  = "integrated.rpca",
                       verbose = FALSE)

seu <- IntegrateLayers(seu, method = CCAIntegration,
                       orig.reduction = "pca",
                       new.reduction  = "integrated.cca",
                       verbose = FALSE)

# Build a temporary UMAP for each reduction (for the comparison figure only).
# Use 1:10 PCs to match paper-replication parameters below.
for (red in c("pca", "harmony", "integrated.rpca", "integrated.cca")) {
  seu <- RunUMAP(seu, reduction = red, dims = 1:10,
                 reduction.name = paste0("umap.", red),
                 verbose = FALSE)
}

# FIGURE 5b - 4-panel integration comparison (colored by sample)
# Build each panel manually so we get clean "UMAP 1 / UMAP 2" axis labels
# (Seurat's auto-named "umappca_1" etc. is unacceptable for publication).
make_int_panel <- function(red_name, title) {
  emb <- Embeddings(seu, red_name)[, 1:2]
  df  <- data.frame(UMAP1 = emb[, 1], UMAP2 = emb[, 2],
                    sample = seu$orig.ident)
  ggplot(df, aes(UMAP1, UMAP2, color = sample)) +
    geom_point(size = 0.25, alpha = 0.7) +
    scale_color_manual(values = SAMPLE_COLORS,
                       labels = c(Day0_Antler = "Velvet",
                                  Day0_Back   = "Back"),
                       name = NULL) +
    labs(title = title, x = "UMAP 1", y = "UMAP 2") +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    theme_umap(8)
}

p_none <- make_int_panel("umap.pca",             "No integration (PCA)")
p_harm <- make_int_panel("umap.harmony",         "Harmony")
p_rpca <- make_int_panel("umap.integrated.rpca", "RPCA")
p_cca  <- make_int_panel("umap.integrated.cca",  "CCA")

fig_int <- (p_none | p_harm) / (p_rpca | p_cca) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Day 0 — integration method comparison",
    subtitle = "UMAP colored by sample (Velvet vs Back). Better integration = uniform color mixing while preserving distinct clusters.",
    theme    = theme(
      plot.title    = element_text(size = 12, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 8, color = "grey30", hjust = 0)
    )
  ) &
  theme(legend.position = "bottom")

save_pub(fig_int, "Integration_comparison", w = 9, h = 8)
# Note: Integration_comparison_split (colored by paper cell type) is generated
# in Section 9, after annotation is done.

saveRDS(seu, file.path(base_dir, "Reindeer_Day0_normalized.rds"))


# 6. Clustering and UMAP
# Paper-replication parameters (Sinha 2022 Day 0): 10 PCs, resolution 0.6.
# Rationale: matches the original publication so UMAP topology is comparable
# to Figure 2A (compact lineage blobs instead of multi-island sub-populations).
# With 10 PCs the fine sub-types (Chondrocyte, Hair_shaft KC, DC, etc.) get
# pulled back into their parent lineage blobs, which is exactly what we want.

INTEGRATION_REDUCTION <- "harmony"
N_PCS                 <- 10        # paper: 10
CLUSTER_RES           <- 0.6       # paper: 0.6

seu <- FindNeighbors(seu, reduction = INTEGRATION_REDUCTION, dims = 1:N_PCS)
seu <- FindClusters(seu, resolution = CLUSTER_RES)
seu <- RunUMAP(seu, reduction = INTEGRATION_REDUCTION, dims = 1:N_PCS,
               reduction.name = "umap")

# FIGURE 6 - UMAP colored by cluster
p_cl <- DimPlot(seu, reduction = "umap", label = TRUE, repel = TRUE,
                label.size = 3, raster = TRUE, pt.size = 0.4) +
  labs(title = sprintf("Clusters (res = %.1f, %d PCs, %s)",
                       CLUSTER_RES, N_PCS, INTEGRATION_REDUCTION),
       x = "UMAP 1", y = "UMAP 2") +
  theme_umap(8) + NoLegend()
save_pub(p_cl, "UMAP_clusters", w = 6, h = 5)

# FIGURE 7 - UMAP colored by sample (post-integration, optimized)
# Optimizations: cell counts in legend, interpretation subtitle,
# equal aspect ratio (coord_fixed), bigger points.
n_velvet <- sum(seu$orig.ident == "Day0_Antler")
n_back   <- sum(seu$orig.ident == "Day0_Back")
emb_smp  <- Embeddings(seu, "umap")
df_smp   <- data.frame(
  UMAP1 = emb_smp[, 1], UMAP2 = emb_smp[, 2],
  tissue = factor(
    ifelse(seu$orig.ident == "Day0_Antler", "Velvet", "Back"),
    levels = c("Velvet", "Back"),
    labels = c(sprintf("Velvet (n = %s)", format(n_velvet, big.mark = ",")),
               sprintf("Back (n = %s)",   format(n_back,   big.mark = ","))))
)
SMP_PAL <- setNames(c(TISSUE_COLORS[["Velvet"]], TISSUE_COLORS[["Back"]]),
                    levels(df_smp$tissue))

p_smp <- ggplot(df_smp, aes(UMAP1, UMAP2, color = tissue)) +
  geom_point(size = 0.5, alpha = 0.8, stroke = 0) +
  scale_color_manual(values = SMP_PAL, name = NULL) +
  coord_fixed() +
  labs(title    = "Sample distribution (post-integration)",
       subtitle = "Velvet and Back cells co-mingle within each cluster - Harmony integration successful",
       x = "UMAP 1", y = "UMAP 2") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme_umap(8) +
  theme(plot.title       = element_text(size = 12, face = "bold", hjust = 0),
        plot.subtitle    = element_text(size = 9, color = "grey30", hjust = 0,
                                        margin = margin(b = 8)),
        legend.position  = c(0.98, 0.98),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = alpha("white", 0.85), color = NA))
save_pub(p_smp, "UMAP_by_sample", w = 7, h = 6)

cat("Number of clusters:", length(unique(seu$seurat_clusters)), "\n")
print(table(seu$seurat_clusters))


# 7. Marker identification
# JoinLayers merges Seurat v5 per-sample layers so FindAllMarkers works.
# Using FindAllMarkers, I identified specific markers for each cluster 
# and annotated 16 cell types by comparing against canonical lineage markers. 

seu <- JoinLayers(seu)
print(Layers(seu))

all_markers <- FindAllMarkers(seu,
                              only.pos        = TRUE,
                              min.pct         = 0.25,
                              logfc.threshold = 0.25)

top10_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) %>%
  ungroup()

write.csv(all_markers,   file.path(base_dir, "Day0_all_markers.csv"),   row.names = FALSE)
write.csv(top10_markers, file.path(base_dir, "Day0_top10_markers.csv"), row.names = FALSE)

saveRDS(seu, file.path(base_dir, "Reindeer_Day0_clustered.rds"))


# 7.5 Schwann-cell diagnostic
# Original paper (Sinha 2022) Day 0 includes a Schwann cluster (SOX10+, MBP+).
# SOX10 is shared with melanocytes, so check MBP/MPZ/PMP22/PLP1 to distinguish.
# If a cluster shows high MBP/MPZ but low MLANA/TYRP1, it is Schwann (not melanocyte).

schwann_markers   <- c("MBP", "MPZ", "PMP22", "PLP1")
melanocyte_markers <- c("MLANA", "TYRP1", "TYR")

# FeaturePlot: Schwann vs Melanocyte distinguishing genes
p_schwann <- FeaturePlot(seu,
                         features = c(schwann_markers, melanocyte_markers, "SOX10"),
                         reduction = "umap", order = TRUE, ncol = 4,
                         raster = TRUE, pt.size = 0.4) &
  scale_color_gradientn(colors = c("grey92", "#5DADE2", "#1A3A66")) &
  labs(x = "UMAP 1", y = "UMAP 2") &
  theme_umap(7) &
  theme(plot.title = element_text(face = "bold.italic", size = 9, hjust = 0.5),
        legend.position = "right",
        legend.key.size = unit(0.3, "cm"))
save_pub(p_schwann, "DIAG_Schwann_vs_Melanocyte", w = 12, h = 6.5)

# Print which clusters are enriched for Schwann markers
schwann_hits <- all_markers[all_markers$gene %in% schwann_markers, ]
cat("\n=== Schwann marker enrichment by cluster ===\n")
if (nrow(schwann_hits) == 0) {
  cat("No Schwann markers found in any cluster's top DEGs.\n")
  cat("Schwann cells may be too few at res=0.5; consider res=0.6 or check manually.\n")
} else {
  print(schwann_hits[, c("cluster", "gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")])
}


# 8. Cell type annotation (two levels: fine + coarse)
# Fine mapping derived from top10 marker review against canonical lineage markers.
# Coarse mapping merges fine types into 9 super-types as requested by PI:
#   Monocyte (Macrophage + DC), Fibroblast (mesenchymal lineages),
#   Keratinocyte (all epithelial including sebaceous), others kept as-is.
# NOTE: cluster numbers may change after integration -- re-check before running.

# --- 8a. Fine cell types (paper-replication run: 10 PCs + res 0.6 -> 14 clusters) ---
# Marker evidence (top10 from FindAllMarkers, Day0_top10_markers.csv):
#   0  COL11A1, DPEP1, FBN2, P4HA3, THBS4         -> Mesenchymal_progenitor
#   1  TCF21, PCOLCE2, MFAP5, A2M, PLAU            -> Fibroblast
#   2  KRT15, LHX2, COL17A1, KRT19, WNT6           -> Basal_keratinocyte (HFSC-flavor)
#   3  KRTDAP, CST6, SBSN, DMKN, KRT6A             -> Suprabasal_keratinocyte
#   4  KRT35, KRT85, KRT89, DLX3                   -> Hair_shaft_keratinocyte
#   5  KERA, MATN4, SFRP4, COL1A1, COL3A1          -> Chondrocyte
#   6  AQP3, FGFBP1, CCL20, HAS3, MMP3             -> Suprabasal_keratinocyte (activated)
#   7  VWF, EMCN, CDH5, PLVAP                      -> Endothelial
#   8  TH, DDC (+ CPA3/TPSB2 mast contamination)   -> Melanocyte (paper-aligned)
#   9  CD3E, CD8A, CTSW, GIMAP7                    -> T_cell
#   10 MKX, COCH, NDNF, IFITM5                     -> Tendon_fibroblast
#   11 CD1B, CD1E, BOLA-DRA/DQB, LYZ               -> Dendritic_cell
#   12 RGS5, MYH11, DES, CNN1                      -> Pericyte_SMC
#   13 C1QA/B/C, FOLR2, FCGR1A                     -> Macrophage
# Schwann (SOX10+/MBP+) not separately resolved at this resolution -- consistent
# with paper Day 0 being Schwann-sparse; flagged in the discussion.
new_ids <- c(
  "0"  = "Mesenchymal_progenitor",
  "1"  = "Fibroblast",
  "2"  = "Basal_keratinocyte",
  "3"  = "Suprabasal_keratinocyte",
  "4"  = "Hair_shaft_keratinocyte",
  "5"  = "Chondrocyte",
  "6"  = "Suprabasal_keratinocyte",
  "7"  = "Endothelial",
  "8"  = "Melanocyte",
  "9"  = "T_cell",
  "10" = "Tendon_fibroblast",
  "11" = "Dendritic_cell",
  "12" = "Pericyte_SMC",
  "13" = "Macrophage"
)
seu <- RenameIdents(seu, new_ids)
seu$celltype_fine <- Idents(seu)
seu$celltype      <- Idents(seu)   # kept for backward-compat with old figures

# --- 8b. Coarse super-types (9 categories, per PI request) ---
coarse_map <- c(
  "Macrophage"              = "Monocyte",
  "Dendritic_cell"          = "Monocyte",
  "Fibroblast"              = "Fibroblast",
  "Chondrocyte"             = "Fibroblast",
  "Tendon_fibroblast"       = "Fibroblast",
  "Mesenchymal_progenitor"  = "Fibroblast",
  "Basal_keratinocyte"      = "Keratinocyte",
  "Suprabasal_keratinocyte" = "Keratinocyte",
  "Hair_shaft_keratinocyte" = "Keratinocyte",
  "Sebaceous_cell"          = "Keratinocyte",
  "Endothelial"             = "Endothelial",
  "Lymphatic_EC"            = "Lymphatic_EC",
  "Pericyte_SMC"            = "Pericyte_SMC",
  "T_cell"                  = "T_cell",
  "Mast_cell"               = "Mast_cell",
  "Melanocyte"              = "Melanocyte"
)
seu$celltype_coarse <- factor(
  unname(coarse_map[as.character(seu$celltype_fine)]),
  levels = unique(coarse_map)
)

# --- 8c. Paper-aligned annotation (Sinha 2022 Cell, Day 0 Fig 2A) ---
# Original Day 0 has 10 labels:
#   Fibroblasts, Basal Keratinocytes, Suprabasal Keratinocyte, Endothelial,
#   Myeloid, Macrophages, T-cell, VSM, Schwann, Melanocyte.
# Mapping rationale:
#   - Basal/Suprabasal kept SEPARATE (paper distinguishes them).
#   - Hair_shaft KC -> Basal (HF lineage is basal-derived).
#   - Sebaceous   -> Suprabasal (terminally differentiated).
#   - Dendritic_cell -> Myeloid; Macrophage -> Macrophages (paper splits these).
#   - Mast_cell      -> Myeloid (myeloid lineage; paper Day 0 has no Mast label).
#   - Lymphatic_EC   -> Endothelial (paper Day 0 does not separate).
#   - Pericyte_SMC   -> VSM (paper terminology).
#   - Schwann is missing from current fine labels -- update new_ids above if
#     the Schwann diagnostic identifies a cluster, then re-run from there.
paper_map <- c(
  "Basal_keratinocyte"      = "Basal Keratinocytes",
  "Hair_shaft_keratinocyte" = "Basal Keratinocytes",
  "Suprabasal_keratinocyte" = "Suprabasal Keratinocyte",
  "Sebaceous_cell"          = "Suprabasal Keratinocyte",
  "Fibroblast"              = "Fibroblasts",
  "Mesenchymal_progenitor"  = "Fibroblasts",
  "Tendon_fibroblast"       = "Fibroblasts",
  "Chondrocyte"             = "Fibroblasts",
  "Endothelial"             = "Endothelial",
  "Lymphatic_EC"            = "Endothelial",
  "Pericyte_SMC"            = "VSM",
  "Macrophage"              = "Macrophages",
  "Dendritic_cell"          = "Myeloid",
  "Mast_cell"               = "Myeloid",
  "T_cell"                  = "T-cell",
  "Melanocyte"              = "Melanocyte",
  "Schwann"                 = "Schwann"   # active only if added to new_ids above
)
seu$celltype_paper <- factor(
  unname(paper_map[as.character(seu$celltype_fine)]),
  levels = unique(paper_map)
)


# 9. Final figures and outputs (both fine and coarse cell-type levels)

# UMAP labeling helper — pastel cloud style:
# slightly larger soft points + translucent alpha gives the "cloud" effect
umap_labeled <- function(obj, group_var, palette, title,
                         label_size = 2.8, legend = FALSE) {
  emb <- Embeddings(obj, "umap")
  df  <- data.frame(UMAP1 = emb[, 1], UMAP2 = emb[, 2],
                    grp   = obj[[group_var, drop = TRUE]])
  centers <- df %>%
    group_by(grp) %>%
    summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2), .groups = "drop")
  p <- ggplot(df, aes(UMAP1, UMAP2, color = grp)) +
    geom_point(size = 0.55, alpha = 0.7, stroke = 0) +
    scale_color_manual(values = palette, name = NULL, drop = FALSE) +
    ggrepel::geom_text_repel(data = centers, aes(label = grp),
                             color = "grey15", size = label_size,
                             bg.color = "white", bg.r = 0.18,
                             fontface = "plain", segment.color = NA,
                             max.overlaps = Inf) +
    labs(title = title, x = "UMAP 1", y = "UMAP 2") +
    guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
    theme_umap(8)
  if (!legend) p <- p + NoLegend()
  p
}

# FIGURE 8a - UMAP with fine cell-type labels
Idents(seu) <- seu$celltype_fine
p_fine <- umap_labeled(seu, "celltype_fine", FINE_COLORS,
                       sprintf("Cell types — fine (%d)",
                               length(unique(seu$celltype_fine))))
save_pub(p_fine, "UMAP_celltype_fine", w = 7.5, h = 6)

# FIGURE 8b - UMAP with coarse super-type labels (9 types, per PI)
Idents(seu) <- seu$celltype_coarse
# Build coarse palette from CELLTYPE_COLORS where possible, fallback to FINE_COLORS
coarse_levels <- levels(seu$celltype_coarse)
coarse_palette <- setNames(
  unname(sapply(coarse_levels, function(x) {
    if (x %in% names(CELLTYPE_COLORS)) CELLTYPE_COLORS[[x]]
    else if (x %in% names(FINE_COLORS)) FINE_COLORS[[x]]
    else "#888888"
  })),
  coarse_levels
)
p_coarse <- umap_labeled(seu, "celltype_coarse", coarse_palette,
                         sprintf("Cell types — coarse (%d)",
                                 length(coarse_levels)))
save_pub(p_coarse, "UMAP_celltype_coarse", w = 7.5, h = 6)

# FIGURE 8b-split - coarse labels, split by sample
emb <- Embeddings(seu, "umap")
df_coarse <- data.frame(UMAP1 = emb[, 1], UMAP2 = emb[, 2],
                        grp = seu$celltype_coarse,
                        tissue = factor(
                          ifelse(seu$orig.ident == "Day0_Antler",
                                 "Velvet", "Back"),
                          levels = c("Back", "Velvet")))
centers_c <- df_coarse %>%
  group_by(tissue, grp) %>%
  summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2), .groups = "drop")
p_coarse_split <- ggplot(df_coarse, aes(UMAP1, UMAP2, color = grp)) +
  geom_point(size = 0.5, alpha = 0.7, stroke = 0) +
  scale_color_manual(values = coarse_palette, name = NULL) +
  ggrepel::geom_text_repel(data = centers_c, aes(label = grp),
                           color = "grey15", size = 2.6, bg.color = "white",
                           bg.r = 0.18, fontface = "plain",
                           segment.color = NA, max.overlaps = Inf) +
  facet_wrap(~ tissue, ncol = 2) +
  labs(title = "Cell types (coarse) — by tissue",
       x = "UMAP 1", y = "UMAP 2") +
  theme_umap(8) + NoLegend() +
  theme(panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.5),
        panel.spacing = unit(1, "lines"))
save_pub(p_coarse_split, "UMAP_celltype_coarse_split", w = 11, h = 5.5)

# FIGURE 8c - UMAP with paper-aligned labels (Sinha 2022 Day 0)
seu$tissue_paper <- factor(
  ifelse(seu$orig.ident == "Day0_Antler", "Velvet", "Back"),
  levels = c("Back", "Velvet")
)
Idents(seu) <- seu$celltype_paper
n_paper <- length(levels(seu$celltype_paper))

paper_palette <- CELLTYPE_COLORS[
  intersect(levels(seu$celltype_paper), names(CELLTYPE_COLORS))]
missing_lab <- setdiff(levels(seu$celltype_paper), names(paper_palette))
if (length(missing_lab) > 0)
  paper_palette <- c(paper_palette, setNames(rep("#888888", length(missing_lab)),
                                             missing_lab))

# Single-panel paper-aligned UMAP
p_paper <- umap_labeled(seu, "celltype_paper", paper_palette,
                        sprintf("Cell types — paper-aligned (%d)", n_paper))
save_pub(p_paper, "UMAP_celltype_paper", w = 7.5, h = 6)

# FIGURE 8c-split - Back | Velvet, colored borders, shared limits (Fig 2A style)
xlims <- range(Embeddings(seu, "umap")[, 1]) + c(-0.5, 0.5)
ylims <- range(Embeddings(seu, "umap")[, 2]) + c(-0.5, 0.5)

make_paper_panel <- function(tissue, border_color) {
  emb <- Embeddings(seu, "umap")
  d   <- data.frame(UMAP1 = emb[, 1], UMAP2 = emb[, 2],
                    grp = seu$celltype_paper,
                    tissue_paper = seu$tissue_paper)
  d   <- d[d$tissue_paper == tissue, ]
  centers <- d %>% group_by(grp) %>%
    summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2), .groups = "drop")
  ggplot(d, aes(UMAP1, UMAP2, color = grp)) +
    geom_point(size = 0.55, alpha = 0.7, stroke = 0) +
    scale_color_manual(values = paper_palette, name = NULL, drop = FALSE) +
    ggrepel::geom_text_repel(data = centers, aes(label = grp),
                             color = "grey15", size = 2.9, bg.color = "white",
                             bg.r = 0.18, fontface = "plain",
                             segment.color = NA, max.overlaps = Inf) +
    coord_cartesian(xlim = xlims, ylim = ylims) +
    labs(title = tissue, x = "UMAP 1", y = "UMAP 2") +
    NoLegend() +
    theme_umap(9) +
    theme(plot.title   = element_text(face = "bold", size = 13, hjust = 0.02,
                                      color = border_color),
          panel.border = element_rect(color = border_color, fill = NA,
                                      linewidth = 1.2))
}

p_back   <- make_paper_panel("Back",   TISSUE_COLORS[["Back"]])
p_velvet <- make_paper_panel("Velvet", TISSUE_COLORS[["Velvet"]])
p_paper_split <- (p_back | p_velvet) +
  plot_annotation(
    title    = "Day 0 cell-type composition — paper-aligned annotation",
    subtitle = sprintf("%d cell types | Harmony integration | %d PCs | res = %.1f",
                       n_paper, N_PCS, CLUSTER_RES),
    theme    = theme(
      plot.title    = element_text(size = 12, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 8.5, color = "grey30", hjust = 0)
    )
  )
save_pub(p_paper_split, "UMAP_celltype_paper_split", w = 11, h = 5.5)

# FIGURE 8d - Integration comparison SPLIT (PI requested to keep)
# 4 integration methods x 2 samples, points colored by paper cell type
# (same palette as main UMAP for cross-figure consistency).
make_split_panel <- function(red_name, title) {
  emb <- Embeddings(seu, red_name)[, 1:2]
  df  <- data.frame(UMAP1   = emb[, 1], UMAP2 = emb[, 2],
                    celltype = seu$celltype_paper,
                    tissue   = factor(
                      ifelse(seu$orig.ident == "Day0_Antler", "Velvet", "Back"),
                      levels = c("Velvet", "Back")))
  ggplot(df, aes(UMAP1, UMAP2, color = celltype)) +
    geom_point(size = 0.55, alpha = 0.9, stroke = 0) +
    scale_color_manual(values = paper_palette, name = NULL, drop = FALSE) +
    facet_wrap(~ tissue, ncol = 2) +
    labs(title = title, x = "UMAP 1", y = "UMAP 2") +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    theme_umap(8) +
    theme(strip.text   = element_text(size = 9, face = "bold"),
          panel.border = element_rect(color = "grey80", fill = NA,
                                      linewidth = 0.3))
}
fig_int_split <- (make_split_panel("umap.pca",             "No integration (PCA)") /
                  make_split_panel("umap.harmony",         "Harmony") /
                  make_split_panel("umap.integrated.rpca", "RPCA") /
                  make_split_panel("umap.integrated.cca",  "CCA")) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Day 0 - Integration comparison (split by sample)",
    subtitle = "Each row is one method. Cell-type structures should look similar across Velvet and Back if integration works.",
    theme = theme(
      plot.title    = element_text(size = 12, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 8, color = "grey30", hjust = 0))) &
  theme(legend.position = "bottom")
save_pub(fig_int_split, "Integration_comparison_split", w = 11, h = 16)

# Reset Idents to fine for the rest of the figures (DotPlot uses fine markers)
Idents(seu) <- seu$celltype_fine

# FIGURE 9 - DotPlot of canonical markers per cell type
key_markers <- c(
  "KRT5","KRT14","KRT15","LHX2",            # Basal KC / HFSC
  "KRTDAP","AQP3","CST6",                   # Suprabasal KC
  "KRT35","KRT85","DLX3",                   # Hair shaft KC
  "COL11A1","SFRP2","DPEP1",                # Mesenchymal progenitor
  "TCF21","PCOLCE2","MFAP5",                # Fibroblast
  "MKX","ASPN","COCH",                      # Tendon fibroblast
  "KERA","MATN4","PODN",                    # Chondrocyte
  "VWF","PLVAP","SELE",                     # Endothelial
  "PROX1","CCL21","MMRN1",                  # Lymphatic EC
  "RGS5","MYH11","DES",                     # Pericyte / SMC
  "MLANA","TYRP1","SOX10",                  # Melanocyte
  "HMGCS2","ACSBG1","SLPI",                 # Sebaceous
  "CD3E","CD8A",                            # T cell
  "C1QA","C1QB","FOLR2",                    # Macrophage
  "CD1B","BOLA-DYA",                        # Dendritic cell
  "CPA3","TPSB2","MS4A2"                    # Mast cell
)
# Use paper-aligned annotation for DotPlot rows (10 types is more readable
# than 13 fine subtypes; this matches the main figure narrative).
Idents(seu) <- seu$celltype_paper
p_dot <- DotPlot(seu, features = unique(key_markers),
                 cols = c("grey90", "#B83668"),   # warm rose, matches Macrophage
                 dot.scale = 4, dot.min = 0.02) +
  labs(x = NULL, y = NULL, title = "Canonical marker expression by cell type") +
  theme_pub(8) +
  theme(
    axis.text.x       = element_text(angle = 45, hjust = 1, size = 7,
                                     face = "italic"),
    axis.text.y       = element_text(size = 8),
    legend.position   = "right",
    legend.box        = "vertical",
    panel.grid.major  = element_line(color = "grey92", linewidth = 0.2)
  ) +
  guides(size  = guide_legend(title = "% expressing",   order = 1),
         color = guide_colorbar(title = "Avg expr (z)", order = 2,
                                barwidth = 0.4, barheight = 4))
save_pub(p_dot, "DotPlot_celltype", w = 13, h = 5.5)

# FIGURE 9b - Top marker heatmap (top 5 per paper-aligned cell type)
# Cell-style: rows = cells (downsampled per type), cols = top genes,
# z-scored expression. Heatmap legend horizontal at top.
Idents(seu) <- seu$celltype_paper
top5_paper <- all_markers %>%
  dplyr::left_join(
    data.frame(cluster = as.character(seq_along(new_ids) - 1),
               fine    = unname(new_ids)),
    by = "cluster") %>%
  dplyr::mutate(paper = unname(paper_map[fine])) %>%
  dplyr::filter(!is.na(paper)) %>%
  dplyr::group_by(paper) %>%
  dplyr::slice_max(order_by = avg_log2FC, n = 5, with_ties = FALSE) %>%
  dplyr::ungroup()

# Downsample to <=80 cells per cell type for legibility
set.seed(42)
keep_cells <- seu@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::group_by(celltype_paper) %>%
  dplyr::slice_sample(n = 80) %>%
  dplyr::pull(cell)

p_hm <- DoHeatmap(seu, cells = keep_cells,
                  features = unique(top5_paper$gene),
                  group.by = "celltype_paper",
                  group.colors = paper_palette[levels(seu$celltype_paper)],
                  size = 2.5, angle = 35, hjust = 0,
                  raster = TRUE) +
  scale_fill_gradientn(colors = c("#2C5985", "white", "#9C2A4D"),
                       na.value = "white",
                       name = "Scaled\nexpression") +
  labs(title = "Top 5 markers per paper-aligned cell type") +
  theme(plot.title  = element_text(size = 11, face = "bold", hjust = 0),
        axis.text.y = element_text(size = 6, face = "italic"),
        legend.position  = "right",
        legend.title     = element_text(size = 8, face = "bold"),
        legend.key.size  = unit(0.35, "cm"))
save_pub(p_hm, "Heatmap_top_markers", w = 12, h = 9)

Idents(seu) <- seu$celltype_fine

# FIGURE 10a - Cell-type proportions per sample (fine)
md <- seu@meta.data %>%
  mutate(tissue = factor(ifelse(orig.ident == "Day0_Antler", "Velvet", "Back"),
                         levels = c("Back", "Velvet")))

p_comp_fine <- ggplot(md, aes(x = tissue, fill = celltype_fine)) +
  geom_bar(position = "fill", color = "white", linewidth = 0.2) +
  scale_fill_manual(values = FINE_COLORS, name = NULL) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Proportion of cells",
       title = "Cell-type composition by tissue — fine") +
  theme_pub(8) +
  theme(legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 7))
save_pub(p_comp_fine, "Composition_by_sample_fine", w = 6, h = 4.5)

# FIGURE 10b - Cell-type proportions per sample (coarse, per PI)
p_comp_coarse <- ggplot(md, aes(x = tissue, fill = celltype_coarse)) +
  geom_bar(position = "fill", color = "white", linewidth = 0.2) +
  scale_fill_manual(values = coarse_palette, name = NULL) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Proportion of cells",
       title = "Cell-type composition by tissue — coarse") +
  theme_pub(8)
save_pub(p_comp_coarse, "Composition_by_sample_coarse", w = 5.5, h = 4.5)

# FIGURE 10c - paper-aligned composition (recommend for paper figure)
p_comp_paper <- ggplot(md, aes(x = tissue, fill = celltype_paper)) +
  geom_bar(position = "fill", color = "white", linewidth = 0.2) +
  scale_fill_manual(values = paper_palette, name = NULL) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Proportion of cells",
       title = "Cell-type composition by tissue — paper-aligned") +
  theme_pub(8)
save_pub(p_comp_paper, "Composition_by_sample_paper", w = 5.5, h = 4.5)

# Cell counts per type x sample (both levels)
cell_counts_fine <- as.data.frame(table(seu$celltype_fine, seu$orig.ident))
colnames(cell_counts_fine) <- c("celltype_fine", "sample", "n_cells")
write.csv(cell_counts_fine,
          file.path(base_dir, "Day0_celltype_fine_counts.csv"),
          row.names = FALSE)

cell_counts_coarse <- as.data.frame(table(seu$celltype_coarse, seu$orig.ident))
colnames(cell_counts_coarse) <- c("celltype_coarse", "sample", "n_cells")
write.csv(cell_counts_coarse,
          file.path(base_dir, "Day0_celltype_coarse_counts.csv"),
          row.names = FALSE)

cell_counts_paper <- as.data.frame(table(seu$celltype_paper, seu$orig.ident))
colnames(cell_counts_paper) <- c("celltype_paper", "sample", "n_cells")
write.csv(cell_counts_paper,
          file.path(base_dir, "Day0_celltype_paper_counts.csv"),
          row.names = FALSE)

cat("\nFine cell types:\n");  print(cell_counts_fine)
cat("\nCoarse cell types:\n"); print(cell_counts_coarse)
cat("\nPaper-aligned cell types:\n"); print(cell_counts_paper)

# Save final annotated object
saveRDS(seu, file.path(base_dir, "Reindeer_Day0_annotated.rds"))


# 10. Fibroblast subclustering (paper Fig 2C reproduction)
# Subset fibroblast lineage, re-run PCA + Harmony + UMAP on the subset.
# UMAP params tuned for a "tree-shaped" connected manifold (paper style):
#   dims = 1:10, min.dist = 0.3, n.neighbors = 200, local.connectivity = 2.

fib_paper_levels <- intersect(c("Fibroblast", "Fibroblasts"),
                              levels(seu$celltype_paper))
fib <- subset(seu, celltype_paper %in% fib_paper_levels)
cat(sprintf("\nFibroblast subset: %d cells\n", ncol(fib)))

# Drop prior dim reductions; re-split RNA layers for re-integration
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
  RunUMAP(reduction = "harmony", dims = 1:10,
          reduction.name = "umap",
          min.dist = 0.3, n.neighbors = 200, spread = 1.0,
          local.connectivity = 2,
          verbose = FALSE)

fib <- JoinLayers(fib)
cat(sprintf("Fibroblast subclusters: %d\n",
            length(unique(fib$seurat_clusters))))

# 12-color deep saturated palette for subclusters
FIB_CLUSTER_COLORS <- c(
  "0"="#C66055", "1"="#AD5818", "2"="#A07C20", "3"="#7C8E2C",
  "4"="#2E7D32", "5"="#2D9686", "6"="#5090C6", "7"="#1F4F94",
  "8"="#9560AE", "9"="#B83668", "10"="#85316C", "11"="#4A4E69"
)

emb_fib   <- Embeddings(fib, "umap")
df_fib_cl <- data.frame(UMAP1 = emb_fib[, 1], UMAP2 = emb_fib[, 2],
                        cluster = fib$seurat_clusters)
centers_fib <- df_fib_cl %>% group_by(cluster) %>%
  summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2), .groups = "drop")

# Paper Fig 2C style: italic "N fibroblasts" tucked above the UMAP arrow
fib_meta <- sprintf("%s fibroblasts", format(ncol(fib), big.mark = ","))
xr <- range(df_fib_cl$UMAP1); yr <- range(df_fib_cl$UMAP2)
meta_x <- xr[1] - diff(xr) * 0.06
meta_y <- yr[1] + diff(yr) * 0.04

p_fib_cl <- ggplot(df_fib_cl, aes(UMAP1, UMAP2, color = cluster)) +
  geom_point(size = 1.1, alpha = 0.95, stroke = 0) +
  scale_color_manual(values = FIB_CLUSTER_COLORS) +
  ggrepel::geom_text_repel(data = centers_fib, aes(label = cluster),
                           color = "black", size = 5, fontface = "bold",
                           bg.color = "white", bg.r = 0.2,
                           segment.color = NA, max.overlaps = Inf) +
  annotate("text", x = meta_x, y = meta_y, label = fib_meta,
           fontface = "italic", size = 3.2, hjust = 0, vjust = 0,
           color = "grey20") +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_umap(9) + NoLegend()
save_pub(p_fib_cl, "UMAP_Fibroblast_subcluster", w = 7, h = 6)

# Save subclustered object
saveRDS(fib, file.path(base_dir, "Reindeer_Day0_Fibroblast_subcluster.rds"))

cat("\n=== Day0 pipeline complete ===\n")
cat("- Cells              :", ncol(seu), "\n")
cat("- Cell types         :", length(unique(seu$celltype)), "\n")
cat("- Fibroblast subclu. :", length(unique(fib$seurat_clusters)),
    sprintf("(%d cells)\n", ncol(fib)))
cat("- Annotated obj      : Reindeer_Day0_annotated.rds\n")
cat("- Fibroblast obj     : Reindeer_Day0_Fibroblast_subcluster.rds\n")
