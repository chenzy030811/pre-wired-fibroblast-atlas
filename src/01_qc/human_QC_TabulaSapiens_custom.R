# ============================================================
# Tabula Sapiens skin QC -- per-donor violins, Before vs After QC
# 2 rows x 3 metrics = 6 panels.
# Same style + same thresholds as Sole-Boldo for direct comparison.
# ============================================================
suppressPackageStartupMessages({
  library(Seurat); library(tidyverse); library(patchwork)
})

base_dir <- "/Users/apple/Downloads/研究生毕业论文/human_data"
raw_dir  <- file.path(base_dir, "Tabula_Sapiens_raw")

# ---- Load 10x-format matrix ----
counts <- ReadMtx(
  mtx      = file.path(raw_dir, "matrix.mtx.gz"),
  cells    = file.path(raw_dir, "barcodes.tsv.gz"),
  features = file.path(raw_dir, "features.tsv.gz"),
  feature.column = 2
)
seu_pre <- CreateSeuratObject(counts, project = "TabulaSapiens",
                              min.cells = 3, min.features = 200)

# Attach donor info from metadata.csv
meta <- read.csv(file.path(raw_dir, "metadata.csv"), row.names = 1)
meta <- meta[colnames(seu_pre), , drop = FALSE]
seu_pre$donor  <- meta$donor_id
seu_pre$tissue <- meta$tissue
seu_pre[["percent.mt"]] <- PercentageFeatureSet(seu_pre, pattern = "^MT-")

# ---- QC thresholds from Tabula Sapiens 2022 paper ----
# Inferred from data ranges in paper-filtered cells:
#   n_genes_by_counts range: 571 - 9987   -> nFeature > 500
#   total_counts range     : 2520 - 114301 -> nCount > 2500
#   pct_counts_mt range    : 0 - 62.86     -> no mt cutoff (use CellBender instead)
TH_FEATURE_LOW  <- 200
TH_FEATURE_HIGH <- 7000
TH_COUNT_LOW    <- 800
TH_COUNT_HIGH   <- 50000
TH_MT           <- 5

seu_post <- subset(seu_pre,
  subset = nFeature_RNA > TH_FEATURE_LOW &
           nFeature_RNA < TH_FEATURE_HIGH &
           nCount_RNA   > TH_COUNT_LOW   &
           nCount_RNA   < TH_COUNT_HIGH  &
           percent.mt   < TH_MT)

cat("Before:", ncol(seu_pre), " After:", ncol(seu_post),
    sprintf(" (%.1f%% retained)\n",
            100 * ncol(seu_post) / ncol(seu_pre)))

# ---- Plot ----
DONOR_COLORS <- c(
  TSP2  = "#E07A5F", TSP10 = "#5BBFA9",
  TSP14 = "#A8C686", TSP21 = "#3D8DAE"
)
donor_levels <- c("TSP2", "TSP10", "TSP14", "TSP21")

theme_pub <- function(base = 9) {
  theme_classic(base_size = base, base_family = "Helvetica") + theme(
    plot.title = element_text(size = base + 2, face = "bold",
                              hjust = 0.5, margin = margin(b = 4)),
    axis.title = element_text(size = base + 1),
    axis.text  = element_text(size = base, color = "black"),
    axis.line  = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.3))
}

qc_violin <- function(df, metric, ylab, low = NULL, high = NULL,
                      y_cap = NULL, title = NULL) {
  p <- ggplot(df, aes(x = donor, y = .data[[metric]], fill = donor)) +
    geom_violin(scale = "width", width = 0.85,
                color = "grey25", linewidth = 0.3, alpha = 0.85) +
    geom_boxplot(width = 0.12, outlier.shape = NA,
                 color = "grey15", fill = "white", linewidth = 0.25) +
    scale_fill_manual(values = DONOR_COLORS, guide = "none") +
    labs(x = NULL, y = ylab, title = title) +
    theme_pub(9) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 0, size = 8))
  # Only draw threshold lines for finite, meaningful values
  if (!is.null(low)  && is.finite(low)  && low  > 0)
    p <- p + geom_hline(yintercept = low,  linetype = "dashed",
                         color = "grey25", linewidth = 0.35)
  if (!is.null(high) && is.finite(high))
    p <- p + geom_hline(yintercept = high, linetype = "dashed",
                         color = "grey25", linewidth = 0.35)
  if (!is.null(y_cap) && is.finite(y_cap))
    p <- p + coord_cartesian(ylim = c(0, y_cap))
  p
}

pre  <- seu_pre@meta.data
post <- seu_post@meta.data
pre$donor  <- factor(pre$donor,  levels = donor_levels)
post$donor <- factor(post$donor, levels = donor_levels)

ycap_f <- quantile(pre$nFeature_RNA, 0.99)
ycap_c <- quantile(pre$nCount_RNA,   0.99)
ycap_m <- as.numeric(quantile(pre$percent.mt, 0.99))   # let data drive cap

p1a <- qc_violin(pre,  "nFeature_RNA", "Genes / cell",
                 TH_FEATURE_LOW, TH_FEATURE_HIGH, ycap_f, "Before QC")
p1b <- qc_violin(pre,  "nCount_RNA",   "UMIs / cell",
                 TH_COUNT_LOW,   TH_COUNT_HIGH,   ycap_c)
p1c <- qc_violin(pre,  "percent.mt",   "Mitochondrial %",
                 NULL, TH_MT, ycap_m)
p2a <- qc_violin(post, "nFeature_RNA", "Genes / cell", title = "After QC")
p2b <- qc_violin(post, "nCount_RNA",   "UMIs / cell")
p2c <- qc_violin(post, "percent.mt",   "Mitochondrial %")

fmt_th <- function(low, high) {
  l <- ifelse(low  == 0   | is.infinite(low),  "-", as.character(low))
  h <- ifelse(is.infinite(high), "-", as.character(high))
  sprintf("%s ~ %s", l, h)
}
hdr <- sprintf(
"Tabula Sapiens 2022 (10x 3' v3, skin only, healthy) | custom thresholds: nFeature %s, nCount %s, %%mt < %s%%
Before QC: %s cells (4 donors)   After QC: %s cells  (%.1f%% retained)",
  fmt_th(TH_FEATURE_LOW, TH_FEATURE_HIGH),
  fmt_th(TH_COUNT_LOW,   TH_COUNT_HIGH),
  ifelse(is.infinite(TH_MT), "-", as.character(TH_MT)),
  format(ncol(seu_pre),  big.mark = ","),
  format(ncol(seu_post), big.mark = ","),
  100 * ncol(seu_post) / ncol(seu_pre))

fig <- (p1a | p1b | p1c) / (p2a | p2b | p2c) +
  plot_annotation(
    title    = "Tabula Sapiens 2022 - QC (uniform 5% mito thresholds for cross-species use)",
    subtitle = hdr,
    theme = theme(
      plot.title    = element_text(size = 12, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 8.5, color = "grey30",
                                   hjust = 0, lineheight = 1.3)))

out_png <- "/Users/apple/Downloads/preview_output/PREVIEW_TabulaSapiens_QC_custom.png"
ggsave(out_png, plot = fig,
       width = 13, height = 6.8, dpi = 300, units = "in",
       bg = "white", device = "png", type = "cairo")

saveRDS(seu_post, file.path(base_dir, "TabulaSapiens_postQC_custom.rds"))
cat("Wrote:", out_png, "\n")
cat("Wrote: human_data/TabulaSapiens_postQC.rds\n")
