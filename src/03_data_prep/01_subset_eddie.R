# ============================================================
# 01 -- Subset Fibroblast + Keratinocyte from 4 species' atlases.
# Runs on Eddie (assumed ~32+ GB RAM available).
# Input : data/*.rds
# Output: subset/{Fibroblast,Keratinocyte}/*.rds  (slim RNA-counts-only)
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
})
options(future.globals.maxSize = 64 * 1024^3)

base <- Sys.getenv("CROSS_ROOT",
                   "/exports/eddie/scratch/YOURUSER/CrossSpecies")
data_dir <- file.path(base, "data")
out_dir  <- file.path(base, "subset")
dir.create(file.path(out_dir, "Fibroblast"),   recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "Keratinocyte"), recursive = TRUE, showWarnings = FALSE)

# Strip everything except RNA counts to keep files small + memory low
slim <- function(seu) {
  if (!"RNA" %in% Assays(seu)) stop("No RNA assay found")
  DefaultAssay(seu) <- "RNA"
  DietSeurat(seu, layers = "counts", assays = "RNA",
             dimreducs = NULL, graphs = NULL)
}

report_save <- function(seu, path, tag) {
  seu <- slim(seu)
  cat(sprintf("  %-16s %6d cells x %5d genes  ->  %s\n",
              tag, ncol(seu), nrow(seu), basename(path)))
  saveRDS(seu, path)
}

# ---------- REINDEER ----------
cat("\n=== Reindeer ===\n")
rdr <- readRDS(file.path(data_dir, "Reindeer_Day0_annotated.rds"))
report_save(subset(rdr, celltype_paper == "Fibroblasts"),
            file.path(out_dir, "Fibroblast",   "Reindeer_Fib.rds"), "Reindeer Fib")
report_save(subset(rdr,
              celltype_paper %in% c("Basal Keratinocytes","Suprabasal Keratinocyte")),
            file.path(out_dir, "Keratinocyte", "Reindeer_KC.rds"),  "Reindeer KC")
rm(rdr); gc()

# ---------- HUMAN (Reynolds) ----------
cat("\n=== Human (Reynolds) ===\n")
hum <- readRDS(file.path(data_dir, "Reynolds_umap.rds"))
report_save(subset(hum, broad_type == "Fibroblast"),
            file.path(out_dir, "Fibroblast",   "Human_Fib.rds"), "Human Fib")
report_save(subset(hum, broad_type == "Keratinocyte"),
            file.path(out_dir, "Keratinocyte", "Human_KC.rds"),  "Human KC")
rm(hum); gc()

# ---------- MUS ----------
cat("\n=== Mus ===\n")
mus <- readRDS(file.path(data_dir, "m.all.final.rds"))
mus_kc_types <- c("IFE suprabasal keratinocytes","IFE basal keratinocytes",
                  "uHF basal keratinocytes",     "uHF suprabasal keratinocytes",
                  "Hair follicle stem cells",     "Proliferating keratinocytes",
                  "Sebaceous glands")
report_save(subset(mus, cell_type == "Dermal fibroblasts"),
            file.path(out_dir, "Fibroblast",   "Mus_Fib.rds"), "Mus Fib")
report_save(subset(mus, cell_type %in% mus_kc_types),
            file.path(out_dir, "Keratinocyte", "Mus_KC.rds"),  "Mus KC")
rm(mus); gc()

# ---------- ACOMYS ----------
cat("\n=== Acomys ===\n")
aco <- readRDS(file.path(data_dir, "a.all.final.rds"))
aco_kc_types <- mus_kc_types  # same annotation vocabulary
report_save(subset(aco, cell_type == "Dermal fibroblasts"),
            file.path(out_dir, "Fibroblast",   "Acomys_Fib.rds"), "Acomys Fib")
report_save(subset(aco, cell_type %in% aco_kc_types),
            file.path(out_dir, "Keratinocyte", "Acomys_KC.rds"),  "Acomys KC")

cat("\n=== DONE. Outputs in", out_dir, "===\n")
