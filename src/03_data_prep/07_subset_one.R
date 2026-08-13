# subset_one.R — process a single species per R process
# Usage:  Rscript subset_one.R  <species>
# where species is one of: reindeer human mus acomys
#
# Each species is processed in its own R session, so memory is fully
# returned to the OS between runs. This lets Mac's 24 GB per-process
# ceiling handle Mus/Acomys without exploding.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1)
  stop("Usage: Rscript subset_one.R  <reindeer|human|mus|acomys>")
species <- tolower(args[1])

suppressPackageStartupMessages({
  library(Seurat)
})
options(future.globals.maxSize = 64 * 1024^3)

base    <- "/Users/apple/Downloads/研究生毕业论文"
out_dir <- file.path(base, "Cross_Species/subset")

slim <- function(seu) {
  DefaultAssay(seu) <- "RNA"
  DietSeurat(seu, layers = "counts", assays = "RNA",
             dimreducs = NULL, graphs = NULL)
}
save_it <- function(seu, path, tag) {
  seu <- slim(seu)
  cat(sprintf("  %-16s %6d cells x %5d genes  ->  %s\n",
              tag, ncol(seu), nrow(seu), basename(path)))
  saveRDS(seu, path)
}

kc_types <- c("IFE suprabasal keratinocytes", "IFE basal keratinocytes",
              "uHF basal keratinocytes",     "uHF suprabasal keratinocytes",
              "Hair follicle stem cells",     "Proliferating keratinocytes",
              "Sebaceous glands")

cat("=== ", toupper(species), " ===\n")

if (species == "reindeer") {
  rdr <- readRDS(file.path(base, "reindeer_data",
                           "Reindeer_Day0_annotated.rds"))
  save_it(subset(rdr, celltype_paper == "Fibroblasts"),
          file.path(out_dir, "Fibroblast",   "Reindeer_Fib.rds"), "Reindeer Fib")
  save_it(subset(rdr, celltype_paper %in%
                       c("Basal Keratinocytes","Suprabasal Keratinocyte")),
          file.path(out_dir, "Keratinocyte", "Reindeer_KC.rds"),  "Reindeer KC")

} else if (species == "human") {
  hum <- readRDS(file.path(base, "human_data", "Reynolds_umap.rds"))
  save_it(subset(hum, broad_type == "Fibroblast"),
          file.path(out_dir, "Fibroblast",   "Human_Fib.rds"), "Human Fib")
  save_it(subset(hum, broad_type == "Keratinocyte"),
          file.path(out_dir, "Keratinocyte", "Human_KC.rds"),  "Human KC")

} else if (species == "mus") {
  mus <- readRDS(file.path(base, "Mus_data", "m.all.final.rds"))
  save_it(subset(mus, cell_type == "Dermal fibroblasts"),
          file.path(out_dir, "Fibroblast",   "Mus_Fib.rds"), "Mus Fib")
  save_it(subset(mus, cell_type %in% kc_types),
          file.path(out_dir, "Keratinocyte", "Mus_KC.rds"),  "Mus KC")

} else if (species == "acomys") {
  aco <- readRDS(file.path(base, "Acomys_data", "a.all.final.rds"))
  save_it(subset(aco, cell_type == "Dermal fibroblasts"),
          file.path(out_dir, "Fibroblast",   "Acomys_Fib.rds"), "Acomys Fib")
  save_it(subset(aco, cell_type %in% kc_types),
          file.path(out_dir, "Keratinocyte", "Acomys_KC.rds"),  "Acomys KC")

} else stop("Unknown species: ", species)

cat("=== ", toupper(species), " DONE ===\n")
