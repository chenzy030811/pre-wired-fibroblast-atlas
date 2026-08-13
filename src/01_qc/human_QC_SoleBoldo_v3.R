# ============================================================
# Sole-Boldo 2020 QC -- per-donor violins, Before vs After QC
# 2 rows x 3 metrics = 6 panels.
# ============================================================
suppressPackageStartupMessages({
  library(Seurat); library(tidyverse); library(patchwork)
})

base_dir <- "/Users/apple/Downloads/研究生毕业论文/human_data"
raw_dir  <- file.path(base_dir, "Sole_Boldo_raw")

# ---- Load raw 10x ----
counts <- ReadMtx(
  mtx      = file.path(raw_dir, "GSE130973_matrix_raw.mtx.gz"),
  cells    = file.path(raw_dir, "GSE130973_barcodes_raw.tsv.gz"),
  features = file.path(raw_dir, "GSE130973_genes_raw.tsv.gz"),
  feature.column = 2
)
seu_pre <- CreateSeuratObject(counts, project = "SoleBoldo",
                              min.cells = 3, min.features = 200)
seu_pre$donor <- sub(".*-", "S", colnames(seu_pre))
seu_pre[["percent.mt"]] <- PercentageFeatureSet(seu_pre, pattern = "^MT-")

# ---- QC thresholds from Sole-Boldo 2020 paper ----
# Extracted from author's saved Seurat object metadata:
#   nCount_RNA range: 800 - 98708  -> lower 800, no upper
#   percent.mito max: 0.04997      -> < 5%
#   nFeature_RNA: no clear threshold (kept down to 40)
TH_FEATURE_LOW  <- 0       # paper used no lower bound
TH_FEATURE_HIGH <- Inf     # paper used no upper bound
TH_COUNT_LOW    <- 800
TH_COUNT_HIGH   <- Inf
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

# 5-sample palette (one color per sequencing lane)
DONOR_COLORS <- c(
  S1 = "#E07A5F", S2 = "#F4A261", S3 = "#A8C686",
  S4 = "#5BBFA9", S5 = "#3D8DAE"
)

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
pre$donor  <- factor(pre$donor,  levels = paste0("S", 1:5))
post$donor <- factor(post$donor, levels = paste0("S", 1:5))

ycap_f <- quantile(pre$nFeature_RNA, 0.99)
ycap_c <- quantile(pre$nCount_RNA,   0.99)
ycap_m <- max(30, TH_MT * 1.6)

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
"Sole-Boldo 2020 (GSE130973) | paper thresholds: nFeature %s, nCount %s, %%mt < %s%%
Before QC: %s cells (5 samples)   After QC: %s cells  (%.1f%% retained)",
  fmt_th(TH_FEATURE_LOW, TH_FEATURE_HIGH),
  fmt_th(TH_COUNT_LOW,   TH_COUNT_HIGH),
  ifelse(is.infinite(TH_MT), "-", as.character(TH_MT)),
  format(ncol(seu_pre),  big.mark = ","),
  format(ncol(seu_post), big.mark = ","),
  100 * ncol(seu_post) / ncol(seu_pre))

fig <- (p1a | p1b | p1c) / (p2a | p2b | p2c) +
  plot_annotation(
    title    = "Sole-Boldo 2020 human skin scRNA-seq - quality control",
    subtitle = hdr,
    theme = theme(
      plot.title    = element_text(size = 12, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 8.5, color = "grey30",
                                   hjust = 0, lineheight = 1.3)))

out_png <- "/Users/apple/Downloads/preview_output/PREVIEW_SoleBoldo_QC.png"
ggsave(out_png, plot = fig,
       width = 13, height = 6.8, dpi = 300, units = "in",
       bg = "white", device = "png", type = "cairo")

saveRDS(seu_post, file.path(base_dir, "SoleBoldo_postQC.rds"))
cat("Wrote:", out_png, "\n")
cat("Wrote: human_data/SoleBoldo_postQC.rds\n")
