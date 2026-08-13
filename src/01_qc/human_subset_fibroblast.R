# ============================================================
# Subset fibroblasts from the full human skin dataset
# Input:  GSE130973_seurat_analysis_lyko 2.rds  (17 clusters, 15,457 cells)
# Output: Human_Fibroblast.rds                 (4 clusters,  ~5,948 cells)
# Fibroblast clusters identified by marker check: 1, 2, 3, 9
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
})

base_dir <- "/Users/apple/Downloads/研究生毕业论文/human_data"

# 1. Load + update
seu <- readRDS(file.path(base_dir, "GSE130973_seurat_analysis_lyko 2.rds"))
seu <- UpdateSeuratObject(seu)
cat("Full dataset:", ncol(seu), "cells,", nrow(seu), "genes\n\n")

# 2. Set cluster column as identity (currently has weird "0_YOUNG" labels)
Idents(seu) <- seu$integrated_snn_res.0.4

# 3. Subset the 4 fibroblast clusters
fib_clusters <- c("1", "2", "3", "9")
fib <- subset(seu, idents = fib_clusters)

cat("Fibroblast subset:", ncol(fib), "cells\n")
cat("Per-cluster breakdown:\n")
print(table(Idents(fib)))
cat("\nPer-age breakdown:\n")
print(table(fib$age))

# 4. Switch back to RNA assay (raw counts) — needed for cross-species integration
DefaultAssay(fib) <- "RNA"

# 5. Save
out_path <- file.path(base_dir, "Human_Fibroblast.rds")
saveRDS(fib, out_path)
cat("\nSaved to:", out_path, "\n")
cat("File size:", round(file.info(out_path)$size / 1024^2, 1), "MB\n")
