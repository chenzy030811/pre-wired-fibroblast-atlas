# ============================================================
# 03a -- For each species subset, rename genes to human symbols
#        using ortholog map, then export 10x-format (mtx + barcodes + features)
#        plus metadata CSV.  Python step (03b) will build the h5ad.
# ============================================================
suppressPackageStartupMessages({
  library(Seurat); library(Matrix)
})

base       <- "/Users/apple/Downloads/研究生毕业论文/Cross_Species"
sub_dir    <- file.path(base, "subset")
out_root   <- file.path(base, "tenx_export")
ortho_path <- file.path(base, "ortholog", "ortholog_map.tsv")

ortho <- read.delim(ortho_path, stringsAsFactors = FALSE)
cat("Ortholog map loaded:", nrow(ortho), "1:1:1:1 genes\n\n")

# species -> ortholog column name
sp_col <- list(Reindeer="reindeer", Human="human", Mus="mouse", Acomys="acomys")

for (comp in c("Fibroblast", "Keratinocyte")) {
  cat("========================================\n")
  cat("Compartment:", comp, "\n")
  cat("========================================\n")
  files <- list.files(file.path(sub_dir, comp), full.names = TRUE)

  for (f in files) {
    tag <- sub("_(Fib|KC)\\.rds$", "", basename(f))  # Reindeer / Human / ...
    species_col <- sp_col[[tag]]
    cat("\n  Loading", basename(f), "( species col =", species_col, ")\n")

    seu <- readRDS(f)
    counts <- GetAssayData(seu, layer = "counts", assay = "RNA")

    # Filter to genes with orthologs, rename to human
    m <- match(rownames(counts), ortho[[species_col]])
    keep <- !is.na(m)
    counts <- counts[keep, , drop = FALSE]
    rownames(counts) <- ortho$human[m[keep]]
    # Drop duplicates (some many:1 possible)
    counts <- counts[!duplicated(rownames(counts)), , drop = FALSE]
    cat("    kept", nrow(counts), "genes (mapped to human symbols)\n")

    out_dir <- file.path(out_root, comp, tag)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    # Write matrix.mtx.gz
    mtx_path <- file.path(out_dir, "matrix.mtx")
    writeMM(counts, mtx_path)
    system(paste("gzip -f", shQuote(mtx_path)))

    # Write barcodes.tsv.gz
    bc_path <- file.path(out_dir, "barcodes.tsv")
    writeLines(colnames(counts), bc_path)
    system(paste("gzip -f", shQuote(bc_path)))

    # Write features.tsv.gz (gene symbol duplicated for id + name columns)
    feat_path <- file.path(out_dir, "features.tsv")
    writeLines(paste(rownames(counts), rownames(counts),
                     "Gene Expression", sep = "\t"), feat_path)
    system(paste("gzip -f", shQuote(feat_path)))

    # Metadata CSV
    meta <- seu@meta.data
    meta$species    <- tag
    meta$compartment <- comp
    write.csv(meta, file.path(out_dir, "metadata.csv"), row.names = TRUE)

    cat("    wrote", out_dir, "\n")
    rm(seu, counts); gc(verbose = FALSE)
  }
}

cat("\n=== DONE ===\n")
