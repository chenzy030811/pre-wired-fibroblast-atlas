# ============================================================
# Reynolds 2021 (E-MTAB-8142) QC -- per-sample violins
# Same 6-panel layout as Sole-Boldo / Tabula Sapiens
# ============================================================
suppressPackageStartupMessages({
  library(Seurat); library(tidyverse); library(patchwork)
})

base_dir <- "/Users/apple/Downloads/研究生毕业论文/human_data"
raw_dir  <- file.path(base_dir, "Reynolds_raw")

# ---- Load sparse MTX (converted from arrayexpress_counts.txt) ----
counts <- ReadMtx(
  mtx      = file.path(raw_dir, "matrix.mtx.gz"),
  cells    = file.path(raw_dir, "barcodes.tsv.gz"),
  features = file.path(raw_dir, "features.tsv.gz"),
  feature.column = 2
)
seu_pre <- CreateSeuratObject(counts, project = "Reynolds",
                              min.cells = 3, min.features = 200)

# ---- Attach metadata ----
meta <- read.delim(file.path(raw_dir, "arrayexpress_metadata.txt"),
                   row.names = 1, stringsAsFactors = FALSE)
meta <- meta[colnames(seu_pre), , drop = FALSE]
seu_pre$sample       <- meta$Sample          # s1, s2, s3, ...
seu_pre$tissue_layer <- meta$Tissue_layer    # Dermis / Epidermis
seu_pre$cell_type    <- meta$Cell_type

seu_pre[["percent.mt"]] <- PercentageFeatureSet(seu_pre, pattern = "^MT-")

cat("Samples found:", paste(unique(seu_pre$sample), collapse = ", "), "\n")
cat("Cells per sample:\n"); print(table(seu_pre$sample))

# ---- QC thresholds from Reynolds 2021 paper ----
# Reverse-engineered from paper-filtered data ranges:
#   nFeature_RNA: 329 - 6502   -> > 200 (data min 329)
#   nCount_RNA  : 1665 - 100102 -> > 1500 (data min 1665)
#   percent.mt  : 0 - 19.98     -> < 20% (cuts exactly at 20)
TH_FEATURE_LOW  <- 200
TH_FEATURE_HIGH <- Inf
TH_COUNT_LOW    <- 1500
TH_COUNT_HIGH   <- Inf
TH_MT           <- 20

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
sample_levels <- sort(unique(seu_pre$sample))
n_samples <- length(sample_levels)
# Cycle through palette if more than 8 samples
palette_base <- c("#E07A5F", "#F4A261", "#E9C46A", "#A8C686",
                  "#5BBFA9", "#3D8DAE", "#5A7AB5", "#9560AE",
                  "#D070B5", "#C66055")
SAMPLE_COLORS <- setNames(
  rep(palette_base, length.out = n_samples), sample_levels)

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
  p <- ggplot(df, aes(x = sample, y = .data[[metric]], fill = sample)) +
    geom_violin(scale = "width", width = 0.85,
                color = "grey25", linewidth = 0.3, alpha = 0.85) +
    geom_boxplot(width = 0.12, outlier.shape = NA,
                 color = "grey15", fill = "white", linewidth = 0.25) +
    scale_fill_manual(values = SAMPLE_COLORS, guide = "none") +
    labs(x = NULL, y = ylab, title = title) +
    theme_pub(9) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 0, size = 8))
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
pre$sample  <- factor(pre$sample,  levels = sample_levels)
post$sample <- factor(post$sample, levels = sample_levels)

ycap_f <- as.numeric(quantile(pre$nFeature_RNA, 0.99))
ycap_c <- as.numeric(quantile(pre$nCount_RNA,   0.99))
ycap_m <- max(30, TH_MT * 6)

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
"Reynolds 2021 (E-MTAB-8142, healthy adult skin) | paper thresholds: nFeature %s, nCount %s, %%mt < %s%%
Before QC: %s cells (%d samples)   After QC: %s cells  (%.1f%% retained)",
  fmt_th(TH_FEATURE_LOW, TH_FEATURE_HIGH),
  fmt_th(TH_COUNT_LOW,   TH_COUNT_HIGH),
  ifelse(is.infinite(TH_MT), "-", as.character(TH_MT)),
  format(ncol(seu_pre),  big.mark = ","),
  n_samples,
  format(ncol(seu_post), big.mark = ","),
  100 * ncol(seu_post) / ncol(seu_pre))

fig <- (p1a | p1b | p1c) / (p2a | p2b | p2c) +
  plot_annotation(
    title    = "Reynolds 2021 human skin scRNA-seq - quality control",
    subtitle = hdr,
    theme = theme(
      plot.title    = element_text(size = 12, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 8.5, color = "grey30",
                                   hjust = 0, lineheight = 1.3)))

out_png <- "/Users/apple/Downloads/preview_output/PREVIEW_Reynolds_QC_paper.png"
ggsave(out_png, plot = fig,
       width = 13, height = 6.8, dpi = 300, units = "in",
       bg = "white", device = "png", type = "cairo")

saveRDS(seu_post, file.path(base_dir, "Reynolds_postQC_paper.rds"))
cat("Wrote:", out_png, "\n")
cat("Wrote: human_data/Reynolds_postQC.rds\n")
