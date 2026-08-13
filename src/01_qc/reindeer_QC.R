# Quality control and Normalization for the Day 0 reindeer dataset
# Data    : GSE168748, Sinha et al. 2022 Cell with 14,540 cell in total
# Samples : Day 0 Antler (velvet, regenerative) + Day 0 Back   (dorsal, scarring)
# Author  : Zhiyi Che n

library(Seurat)
library(tidyverse)
library(patchwork)

# 1. Paths, samples, thresholds (single source of truth) 
base_dir <- "/Users/apple/Downloads/研究生毕业论文/reindeer_data"
samples  <- c("Day0_Antler", "Day0_Back")

# I added lower bounds to prevent low-feature cells from forming an artefact cluster downstream. 
TH_FEATURE_LOW  <- 800
TH_FEATURE_HIGH <- 4500
TH_COUNT_LOW    <- 1500
TH_COUNT_HIGH   <- 15000
TH_MT           <- 5

# Reindeer is annotated against the bovine reference, so mitochondrial
# genes use gene-symbol naming (ND/COX/ATP/CYTB) rather than the usual
# "MT-" prefix. This regex covers all expected naming variants.
MT_PATTERN <- "^MT-|^mt-|^Mt-|^ND[1-6]$|^ND4L$|^COX[1-3]$|^ATP[68]$|^CYTB$"

SAMPLE_COLORS <- c(Day0_Antler = "#E07A5F",
                   Day0_Back   = "#3D8DAE")


# 2. Load each sample, light pre-filter, merge 
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


# 3. Apply QC filter 
seu_post <- subset(seu_pre,
                   subset = nFeature_RNA > TH_FEATURE_LOW &
                     nFeature_RNA < TH_FEATURE_HIGH &
                     nCount_RNA   > TH_COUNT_LOW &
                     nCount_RNA   < TH_COUNT_HIGH &
                     percent.mt   < TH_MT)

cat("Antler:", sum(seu_post$orig.ident == "Day0_Antler"), "\n")
cat("Back:  ", sum(seu_post$orig.ident == "Day0_Back"),   "\n")
cat("Total: ", ncol(seu_post), "/", ncol(seu_pre),
    sprintf("(%.1f%%)\n", 100 * ncol(seu_post) / ncol(seu_pre)))

# Save the cleaned object for downstream analysis
saveRDS(seu_post, file.path(base_dir, "Reindeer_Day0_postQC.rds"))


# FIGURE 1 — QC report (before/after violins + scatters)
# Helper: 3-panel violin (nFeature, nCount, percent.mt) with threshold lines
make_vln_annotated <- function(obj, title) {
  p_feat <- VlnPlot(obj, "nFeature_RNA", group.by = "orig.ident", pt.size = 0) +
    geom_hline(yintercept = c(TH_FEATURE_LOW, TH_FEATURE_HIGH),
               linetype = "dashed", color = "red") + NoLegend()
  p_cnt  <- VlnPlot(obj, "nCount_RNA", group.by = "orig.ident", pt.size = 0) +
    geom_hline(yintercept = c(TH_COUNT_LOW, TH_COUNT_HIGH),
               linetype = "dashed", color = "red") + NoLegend()
  p_mt   <- VlnPlot(obj, "percent.mt", group.by = "orig.ident", pt.size = 0) +
    geom_hline(yintercept = TH_MT, linetype = "dashed", color = "red") + NoLegend()
  (p_feat | p_cnt | p_mt) + plot_annotation(title = title)
}

# Helper: 2-panel scatter with threshold lines on both axes
make_scatter <- function(obj, title) {
  p1 <- FeatureScatter(obj, "nCount_RNA", "percent.mt",
                       group.by = "orig.ident") +
    geom_hline(yintercept = TH_MT, linetype = "dashed", color = "red") +
    geom_vline(xintercept = c(TH_COUNT_LOW, TH_COUNT_HIGH),
               linetype = "dashed", color = "red") +
    ggtitle("UMI vs. percent.mt")
  p2 <- FeatureScatter(obj, "nCount_RNA", "nFeature_RNA",
                       group.by = "orig.ident") +
    geom_hline(yintercept = c(TH_FEATURE_LOW, TH_FEATURE_HIGH),
               linetype = "dashed", color = "red") +
    geom_vline(xintercept = c(TH_COUNT_LOW, TH_COUNT_HIGH),
               linetype = "dashed", color = "red") +
    ggtitle("UMI vs. nFeature")
  (p1 | p2) + plot_annotation(title = title)
}

# Summary text block (all numbers are computed dynamically)
pct_removed <- 100 * (1 - ncol(seu_post) / ncol(seu_pre))
summary_text <- sprintf(
  "QC thresholds applied:
  nFeature_RNA  in  (%d, %d)
  nCount_RNA    in  (%d, %d)
  percent.mt    <   %d%%

Before QC:  %d cells   (Antler %d  |  Back %d)
After  QC:  %d cells   (Antler %d  |  Back %d)
Removed:    %d cells   (%.1f%%)",
  TH_FEATURE_LOW, TH_FEATURE_HIGH,
  TH_COUNT_LOW,   TH_COUNT_HIGH,
  TH_MT,
  ncol(seu_pre),  sum(seu_pre$orig.ident  == "Day0_Antler"),
  sum(seu_pre$orig.ident  == "Day0_Back"),
  ncol(seu_post), sum(seu_post$orig.ident == "Day0_Antler"),
  sum(seu_post$orig.ident == "Day0_Back"),
  ncol(seu_pre) - ncol(seu_post), pct_removed)

p_text <- ggplot() +
  annotate("text", x = 0, y = 1, hjust = 0, vjust = 1,
           label = summary_text, family = "mono", size = 4) +
  xlim(0, 1) + ylim(0, 1) + theme_void()

final_fig <- p_text /
  make_vln_annotated(seu_pre,  "Before QC") /
  make_vln_annotated(seu_post, "After QC")  /
  make_scatter(seu_pre,  "Before QC — scatter") /
  make_scatter(seu_post, "After QC — scatter") +
  plot_layout(heights = c(1, 2, 2, 2, 2))

ggsave(file.path(base_dir, "QC_report_full.png"),
       plot = final_fig, width = 13, height = 18, dpi = 300)

# FIGURE 2 — Per-sample violins with jittered points
# Saved both as individual PNGs and as a combined 2x3 panel
# Helper: one publication-style violin (background grey = filtered,
# colored overlay = kept, jittered points on top, threshold lines)
make_qc_violin <- function(obj, sample_name, metric,
                           low_th = NULL, high_th = NULL, y_label) {
  
  df <- FetchData(obj, c("orig.ident", metric)) %>%
    filter(orig.ident == sample_name)
  color <- SAMPLE_COLORS[[sample_name]]
  
  # Per-cell filter status under this metric's thresholds
  kept_mask <- rep(TRUE, nrow(df))
  if (!is.null(low_th))  kept_mask <- kept_mask & df[[metric]] >= low_th
  if (!is.null(high_th)) kept_mask <- kept_mask & df[[metric]] <= high_th
  df$status <- ifelse(kept_mask, "kept", "filtered")
  df_kept   <- df[kept_mask, ]
  
  # Trim y-axis so long tails don't stretch the plot
  y_max <- if (!is.null(high_th)) high_th * 1.1
  else as.numeric(quantile(df[[metric]], 0.995, na.rm = TRUE))
  
  pct_kept <- round(100 * nrow(df_kept) / nrow(df), 1)
  
  ggplot(df, aes(x = orig.ident, y = .data[[metric]])) +
    geom_violin(fill = "grey80", color = "grey55",
                scale = "width", linewidth = 0.3, width = 0.85) +
    geom_violin(data = df_kept,
                fill = color, color = "black",
                scale = "width", linewidth = 0.4, alpha = 0.85, width = 0.85) +
    geom_jitter(aes(color = status),
                width = 0.32, size = 0.25, alpha = 0.45, stroke = 0) +
    scale_color_manual(values = c(kept = "black", filtered = "grey50")) +
    { if (!is.null(low_th))
      geom_hline(yintercept = low_th,  linetype = "dashed",
                 color = "grey25", linewidth = 0.4) } +
    { if (!is.null(high_th))
      geom_hline(yintercept = high_th, linetype = "dashed",
                 color = "grey25", linewidth = 0.4) } +
    coord_cartesian(ylim = c(0, y_max)) +
    labs(title    = sprintf("%s — %s", sample_name, metric),
         subtitle = sprintf("kept %d / %d cells  (%.1f%%)",
                            nrow(df_kept), nrow(df), pct_kept),
         x = NULL, y = y_label) +
    theme_classic(base_size = 11, base_family = "Helvetica") +
    theme(
      plot.title        = element_text(size = 12, face = "bold", hjust = 0.5,
                                       margin = margin(0, 0, 0, 0)),
      plot.subtitle     = element_text(size = 9, color = "grey35", hjust = 0.5,
                                       margin = margin(0, 0, 0, 0)),
      axis.title.y      = element_text(size = 11, margin = margin(r = 6)),
      axis.text         = element_text(size = 10, color = "black"),
      axis.ticks        = element_line(linewidth = 0.3, color = "black"),
      axis.ticks.length = unit(-0.12, "cm"),
      axis.text.x       = element_blank(),
      axis.ticks.x      = element_blank(),
      axis.line         = element_line(linewidth = 0.4, color = "black"),
      legend.position   = "none",
      plot.margin       = margin(6, 10, 6, 6)
    )
}

# Metric specs (NULL = no threshold on that side)
metric_specs <- list(
  list(name = "nFeature_RNA", low = TH_FEATURE_LOW, high = TH_FEATURE_HIGH, ylab = "Genes / cell"),
  list(name = "nCount_RNA",   low = TH_COUNT_LOW,   high = TH_COUNT_HIGH,   ylab = "UMIs / cell"),
  list(name = "percent.mt",   low = NULL,           high = TH_MT,           ylab = "Mitochondrial %")
)

# Build all 6 plots once, save each individually and combine
plots <- list()
for (s in samples) {
  for (m in metric_specs) {
    p <- make_qc_violin(seu_pre, s, m$name, m$low, m$high, m$ylab)
    plots[[length(plots) + 1]] <- p
    
    ggsave(file.path(base_dir, sprintf("VlnPlot_%s_%s.png", s, m$name)),
           plot = p, width = 3, height = 4.5, dpi = 300)
  }
}

# Combined 2 rows (samples) x 3 cols (metrics)
combined <- wrap_plots(plots, ncol = 3) +
  plot_annotation(
    title    = "Reindeer Day 0 — Per-sample QC violins",
    subtitle = sprintf(
      "Final after all filters:  Antler %d  |  Back %d  (total %d / %d, %.1f%% retained)",
      sum(seu_post$orig.ident == "Day0_Antler"),
      sum(seu_post$orig.ident == "Day0_Back"),
      ncol(seu_post), ncol(seu_pre),
      100 * ncol(seu_post) / ncol(seu_pre)),
    theme = theme(
      plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, color = "grey35", hjust = 0.5))
  )

ggsave(file.path(base_dir, "QC_violins_combined.png"),
       plot = combined, width = 12, height = 9, dpi = 300)

cat("\nAll figures and Reindeer_Day0_postQC.rds saved to:", base_dir, "\n")


# After filtering, 9,412 out of 14,540 cells (64.7%) were retained: 3,175 from Antler and 6237 from Back. 
# From the violin plots, both nCount and percent.mt show a shaped distribution with a central peak and very few at the extremes;
# nFeature exhibits a slight shoulder on the lower end due to a dense cluster of granulocytes, 
# which is an inherent limitation of the data.

