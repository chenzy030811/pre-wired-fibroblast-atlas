# Reindeer Day0 (Antler + Back) — Normalization & HVG selection
# After QC, data were log-normalized (scale factor = 10,000) to correct for sequencing depth differences. 
# High variable genes (HVGs) were selected using vst with a threshold of 2,000, and then z-score scaled in preparation for PCA.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

# paths
base_dir <- "/Users/apple/Downloads/研究生毕业论文/reindeer_data"
in_rds   <- file.path(base_dir, "Reindeer_Day0_postQC.rds")
out_rds  <- file.path(base_dir, "Reindeer_Day0_normalized.rds")

# load post-QC object
seu <- readRDS(in_rds)
cat("Loaded object:", ncol(seu), "cells x", nrow(seu), "genes\n")
cat("Samples:\n"); print(table(seu$orig.ident))

# Normalize: log( 10000 * count / total + 1 ) 
seu <- NormalizeData(seu,
                     normalization.method = "LogNormalize",
                     scale.factor         = 1e4)

# Highly variable genes (top 2000, vst) 
seu <- FindVariableFeatures(seu,
                            selection.method = "vst",
                            nfeatures        = 2000)

top10 <- head(VariableFeatures(seu), 10)
cat("\nTop 10 HVGs:\n"); print(top10)

# HVG mean–variance plot (avoid -Inf from zero-variance genes)
p1 <- VariableFeaturePlot(seu) 
p2 <- LabelPoints(plot   = p1,
                  points = top10,
                  repel  = TRUE,
                  xnudge = 0,
                  ynudge = 0)
ggsave(file.path(base_dir, "HVG_plot.png"),
       plot = p2, width = 8, height = 5, dpi = 300)

# Scale data 
# features = rownames(seu): scale ALL genes so DoHeatmap/DotPlot
#   on non-HVG markers later won't fail.
# vars.to.regress: remove technical effects of library size & mito%.
seu <- ScaleData(seu,
                 features        = rownames(seu),
                 vars.to.regress = c("nCount_RNA", "percent.mt"))

# Save 
saveRDS(seu, out_rds)

# Summary
cat("\n=== Normalization complete ===\n")
cat("- Cells     :", ncol(seu), "\n")
cat("- Genes     :", nrow(seu), "\n")
cat("- HVGs      :", length(VariableFeatures(seu)), "\n")
cat("- Top 10 HVG:", paste(top10, collapse = ", "), "\n")
cat("- Saved to  :", out_rds, "\n")

while (!is.null(dev.list())) dev.off()

top10 <- head(VariableFeatures(seu), 10)
p1 <- VariableFeaturePlot(seu)
p2 <- LabelPoints(p1, points = top10, repel = TRUE,
                  xnudge = 0, ynudge = 0)
suppressWarnings(
  ggsave(file.path(base_dir, "HVG_plot.png"), p2,
         width = 8, height = 5, dpi = 300)
)

seu <- RunPCA(seu, npcs = 30)
ElbowPlot(seu, ndims = 30)
png(file.path(base_dir, "PCA_by_sample.png"),
    width = 7, height = 5, units = "in", res = 300)
print(DimPlot(seu, reduction = "pca", group.by = "orig.ident"))
dev.off()

saveRDS(seu, file.path(base_dir, "Reindeer_Day0_normalized.rds"))

# Initial PCA shows clean biological signals (epithelial, fibroblast, immune, endothelial, mast cells) 
# with minimal batch effect between the two samples, so no integration is needed. Next step: clustering and UMAP visualization.


