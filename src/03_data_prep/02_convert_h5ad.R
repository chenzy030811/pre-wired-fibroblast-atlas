# ============================================================
# 03 -- Convert each species subset -> h5ad, rename genes to human,
#       keep only the 1:1:1:1 ortholog gene set, concatenate.
# Output: h5ad/{Fibroblast,Keratinocyte}/{species}.h5ad
#         + merged_{Fibroblast,Keratinocyte}.h5ad  (input for scVI)
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  if (!requireNamespace("MuDataSeurat", quietly = TRUE)) {
    if (!requireNamespace("remotes", quietly = TRUE))
      install.packages("remotes", repos = "https://cloud.r-project.org")
    remotes::install_github("PMBio/MuDataSeurat", quiet = TRUE)
  }
  library(MuDataSeurat)
})

base <- Sys.getenv("CROSS_ROOT",
                   "/exports/eddie/scratch/YOURUSER/CrossSpecies")
sub_dir  <- file.path(base, "subset")
ortho    <- read.delim(file.path(base, "ortholog", "ortholog_map.tsv"),
                        stringsAsFactors = FALSE)
h5ad_dir <- file.path(base, "h5ad")

for (comp in c("Fibroblast", "Keratinocyte")) {
  cat("\n============================================================\n")
  cat("Compartment:", comp, "\n")
  cat("============================================================\n")

  dir.create(file.path(h5ad_dir, comp), recursive = TRUE, showWarnings = FALSE)
  files <- list.files(file.path(sub_dir, comp), full.names = TRUE)

  merged_list <- list()

  for (f in files) {
    tag <- sub("_(Fib|KC)\\.rds$", "", basename(f))   # Reindeer / Human / Mus / Acomys
    species_col <- tolower(tag)
    if (species_col == "reindeer") species_col <- "reindeer"

    cat("\n  Loading:", basename(f), "\n")
    seu <- readRDS(f)
    counts <- GetAssayData(seu, layer = "counts", assay = "RNA")

    # Ortholog rename: filter counts to only genes with a 1:1 human ortholog
    keep <- rownames(counts) %in% ortho[[species_col]]
    counts <- counts[keep, , drop = FALSE]
    idx <- match(rownames(counts), ortho[[species_col]])
    rownames(counts) <- ortho$human[idx]
    counts <- counts[!is.na(rownames(counts)), , drop = FALSE]
    counts <- counts[!duplicated(rownames(counts)), , drop = FALSE]
    cat("    ->", nrow(counts), "genes kept (mapped to human symbols)\n")

    # Add species column to metadata
    meta <- seu@meta.data
    meta$species <- tag

    # Build a slim Seurat with the renamed matrix
    seu2 <- CreateSeuratObject(counts = counts, meta.data = meta)
    seu2$species <- tag
    seu2$compartment <- comp

    # Save per-species h5ad
    out_h5 <- file.path(h5ad_dir, comp, paste0(tag, ".h5ad"))
    MuDataSeurat::WriteH5AD(seu2, file = out_h5, assay = "RNA")
    cat("    wrote:", out_h5, "\n")

    merged_list[[tag]] <- seu2
    rm(seu, seu2, counts); gc()
  }

  # Merge across species (union of cells, intersection of ortholog genes)
  cat("\n  Merging 4 species...\n")
  merged <- merge(merged_list[[1]], y = merged_list[-1],
                  add.cell.ids = names(merged_list))
  merged <- JoinLayers(merged)
  cat("  Merged dims:", nrow(merged), "genes x", ncol(merged), "cells\n")

  out_merged <- file.path(h5ad_dir, paste0("merged_", comp, ".h5ad"))
  MuDataSeurat::WriteH5AD(merged, file = out_merged, assay = "RNA")
  cat("  wrote:", out_merged, "\n")
}

cat("\n=== DONE ===\n")
