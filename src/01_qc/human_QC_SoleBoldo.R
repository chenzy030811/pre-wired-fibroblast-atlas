# ============================================================
# QC of Sole-Boldo 2020 (GSE130973) from RAW 10x matrices
# Goal: produce QC distributions to compare with Reynolds 2021
# Pipeline mirrors the reindeer QC step (same logic, same thresholds)
# ============================================================
suppressPackageStartupMessages({
  library(Seurat); library(tidyverse); library(patchwork)
})

# Manual skewness / kurtosis (avoid moments package dependency)
skewness <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x); m <- mean(x); s <- sd(x)
  sum(((x - m) / s)^3) / n
}
kurtosis <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x); m <- mean(x); s <- sd(x)
  sum(((x - m) / s)^4) / n
}

base_dir <- "/Users/apple/Downloads/研究生毕业论文/human_data"
raw_dir  <- file.path(base_dir, "Sole_Boldo_raw")

# ---- 1. Load raw 10x matrix (custom file names) ----
cat("Loading raw 10x matrix...\n")
counts <- ReadMtx(
  mtx      = file.path(raw_dir, "GSE130973_matrix_raw.mtx.gz"),
  cells    = file.path(raw_dir, "GSE130973_barcodes_raw.tsv.gz"),
  features = file.path(raw_dir, "GSE130973_genes_raw.tsv.gz"),
  feature.column = 2   # gene symbol column
)
cat("  Raw matrix size :", nrow(counts), "genes x", ncol(counts), "droplets\n")

# ---- 2. Light filter (remove obvious empty droplets) ----
# NOT a QC step; this just trims empty drops so distributions aren't dominated
# by empty water droplets. Same as the reindeer setup.
seu <- CreateSeuratObject(counts,
                          project = "SoleBoldo",
                          min.cells    = 3,
                          min.features = 200)
cat("  After empty-drop trim:", ncol(seu), "cells\n")

# ---- 3. Parse donor identity from barcode suffix (-1, -2, ..., -10) ----
seu$donor <- sub(".*-", "S", colnames(seu))   # "AAACC-3" -> "S3"
cat("\nCells per donor:\n"); print(table(seu$donor))

# ---- 4. Compute QC metrics ----
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

cat("\nQC summary (median across all cells):\n")
qc_summary <- seu@meta.data %>%
  summarise(
    median_nFeature = median(nFeature_RNA),
    median_nCount   = median(nCount_RNA),
    median_percent_mt = median(percent.mt),
    cells = n()
  )
print(qc_summary)

# ---- 5. Normality tests + shape stats per metric ----
metrics <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
shape_stats <- lapply(metrics, function(m) {
  x <- seu@meta.data[[m]]
  x <- x[!is.na(x)]
  log_x <- log1p(x)             # log-transform handles long-tail counts
  ks    <- ks.test(scale(log_x)[, 1], "pnorm")
  data.frame(
    metric         = m,
    skewness_raw   = round(skewness(x), 2),
    kurtosis_raw   = round(kurtosis(x), 2),
    skewness_log   = round(skewness(log_x), 2),
    kurtosis_log   = round(kurtosis(log_x), 2),
    ks_pvalue_log  = signif(ks$p.value, 3)
  )
}) %>% bind_rows()

cat("\nShape statistics (closer to skew=0, kurt=3 = more normal):\n")
print(shape_stats)
write.csv(shape_stats,
          file.path(base_dir, "SoleBoldo_QC_shape_stats.csv"),
          row.names = FALSE)

# ---- 6. QC plots (raw + log scale + density overlay) ----
SAMPLE_COL <- "#E07A5F"          # consistent with reindeer Velvet color

make_qc <- function(df, metric, label, log_scale = FALSE) {
  p <- ggplot(df, aes(x = .data[[metric]])) +
    geom_density(fill = SAMPLE_COL, color = "grey20", alpha = 0.7,
                 linewidth = 0.3) +
    labs(x = label, y = "Density",
         title = paste0(label, ifelse(log_scale, "  (log10)", ""))) +
    theme_classic(base_size = 9, base_family = "Helvetica") +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5))
  if (log_scale) p <- p + scale_x_log10()
  p
}

md <- seu@meta.data

p1a <- make_qc(md, "nFeature_RNA", "Genes / cell")
p1b <- make_qc(md, "nFeature_RNA", "Genes / cell", log_scale = TRUE)
p2a <- make_qc(md, "nCount_RNA",   "UMIs / cell")
p2b <- make_qc(md, "nCount_RNA",   "UMIs / cell", log_scale = TRUE)
p3a <- make_qc(md, "percent.mt",   "Mitochondrial %")
p3b <- ggplot(md, aes(x = log1p(percent.mt))) +
  geom_density(fill = SAMPLE_COL, color = "grey20", alpha = 0.7) +
  labs(x = "log(1 + %mt)", y = "Density",
       title = "Mitochondrial %  (log1p)") +
  theme_classic(base_size = 9, base_family = "Helvetica") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5))

# Per-donor violin (10 donors)
md$donor <- factor(md$donor, levels = paste0("S", 1:10))
pv1 <- ggplot(md, aes(donor, nFeature_RNA, fill = donor)) +
  geom_violin(scale = "width", color = "grey20", alpha = 0.8, linewidth = 0.3) +
  geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA, linewidth = 0.3) +
  labs(x = NULL, y = "Genes / cell", title = "Per-donor: nFeature_RNA") +
  theme_classic(base_size = 9) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, size = 8),
        plot.title = element_text(face = "bold"))

pv2 <- ggplot(md, aes(donor, nCount_RNA, fill = donor)) +
  geom_violin(scale = "width", color = "grey20", alpha = 0.8, linewidth = 0.3) +
  geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA, linewidth = 0.3) +
  labs(x = NULL, y = "UMIs / cell", title = "Per-donor: nCount_RNA") +
  theme_classic(base_size = 9) + NoLegend() +
  theme(plot.title = element_text(face = "bold"))

pv3 <- ggplot(md, aes(donor, percent.mt, fill = donor)) +
  geom_violin(scale = "width", color = "grey20", alpha = 0.8, linewidth = 0.3) +
  geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA, linewidth = 0.3) +
  labs(x = NULL, y = "Mitochondrial %", title = "Per-donor: percent.mt") +
  coord_cartesian(ylim = c(0, 20)) +
  theme_classic(base_size = 9) + NoLegend() +
  theme(plot.title = element_text(face = "bold"))

stat_text <- sprintf(
  "Sole-Boldo (GSE130973) - raw, %s cells, %d donors\nMedian nFeature %s, nCount %s, %%mt %.2f%%\nSkewness (log-scale) - nFeature %.2f, nCount %.2f, %%mt %.2f",
  format(ncol(seu), big.mark = ","), 10,
  format(qc_summary$median_nFeature, big.mark = ","),
  format(qc_summary$median_nCount, big.mark = ","),
  qc_summary$median_percent_mt,
  shape_stats$skewness_log[shape_stats$metric == "nFeature_RNA"],
  shape_stats$skewness_log[shape_stats$metric == "nCount_RNA"],
  shape_stats$skewness_log[shape_stats$metric == "percent.mt"]
)

fig <- (p1a | p2a | p3a) / (p1b | p2b | p3b) / (pv1 | pv2 | pv3) +
  plot_annotation(
    title    = "Sole-Boldo 2020 (GSE130973) - QC distribution from raw 10x",
    subtitle = stat_text,
    theme = theme(
      plot.title    = element_text(size = 12, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 8.5, color = "grey30", hjust = 0,
                                   lineheight = 1.3)))

ggsave("/tmp/PREVIEW_SoleBoldo_QC.png", plot = fig,
       width = 12, height = 10, dpi = 300, units = "in",
       bg = "white", device = "png", type = "cairo")

# ---- 7. Save QC'd Seurat object for later cross-dataset comparison ----
saveRDS(seu, file.path(base_dir, "SoleBoldo_raw_qc.rds"))

cat("\nWrote:\n")
cat("  /tmp/PREVIEW_SoleBoldo_QC.png\n")
cat("  human_data/SoleBoldo_QC_shape_stats.csv\n")
cat("  human_data/SoleBoldo_raw_qc.rds\n")
cat("\n=== DONE ===\n")
