# ============================================================
# Step 1 -- Subset Fibroblast + Keratinocyte from 4 species
# Output: 8 RDS files under Cross_Species/{Fibroblast,Keratinocyte}/
# ============================================================
Sys.setenv(R_MAX_VSIZE = 64e9)   # 64 GB ceiling (macOS respects this)
suppressPackageStartupMessages({
  library(Seurat); library(dplyr)
})
options(future.globals.maxSize = 32 * 1024^3)

# Strip heavy layers before saving — we only need RNA counts for downstream
# ortholog + scVI. SCT scale.data is what blows up memory.
slim_object <- function(seu) {
  DefaultAssay(seu) <- "RNA"           # switch away from SCT before pruning
  seu <- DietSeurat(seu,
                    layers    = "counts",
                    assays    = "RNA",
                    dimreducs = NULL,
                    graphs    = NULL)
  seu
}

base_dir <- "/Users/apple/Downloads/研究生毕业论文"
out_dir  <- file.path(base_dir, "Cross_Species")
dir.create(file.path(out_dir, "Fibroblast"),   recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "Keratinocyte"), recursive = TRUE, showWarnings = FALSE)

save_subset <- function(seu, out_path, species) {
  seu <- slim_object(seu)
  cat(sprintf("  %s : %d cells x %d genes  ->  %s\n",
              species, ncol(seu), nrow(seu), basename(out_path)))
  saveRDS(seu, out_path)
}

# =====================================================
# REINDEER
# =====================================================
cat("\n=== Reindeer ===\n")
rdr <- readRDS(file.path(base_dir, "reindeer_data",
                         "Reindeer_Day0_annotated.rds"))

# Fibroblast (paper label; already includes MP, Chondrocyte, Tendon_fib)
rdr_fib <- subset(rdr, celltype_paper == "Fibroblasts")
save_subset(rdr_fib,
            file.path(out_dir, "Fibroblast", "Reindeer_Fib.rds"),
            "Reindeer Fib")

# Keratinocyte (paper label; both Basal + Suprabasal)
rdr_kc <- subset(rdr,
  celltype_paper %in% c("Basal Keratinocytes", "Suprabasal Keratinocyte"))
save_subset(rdr_kc,
            file.path(out_dir, "Keratinocyte", "Reindeer_KC.rds"),
            "Reindeer KC")

# =====================================================
# HUMAN (Reynolds)
# =====================================================
cat("\n=== Human (Reynolds) ===\n")
hum <- readRDS(file.path(base_dir, "human_data", "Reynolds_umap.rds"))

hum_fib <- subset(hum, broad_type == "Fibroblast")
save_subset(hum_fib,
            file.path(out_dir, "Fibroblast", "Human_Fib.rds"),
            "Human Fib")

hum_kc <- subset(hum, broad_type == "Keratinocyte")
save_subset(hum_kc,
            file.path(out_dir, "Keratinocyte", "Human_KC.rds"),
            "Human KC")

# =====================================================
# MUS
# =====================================================
cat("\n=== Mus ===\n")
mus <- readRDS(file.path(base_dir, "Mus_data", "m.all.final.rds"))

mus_fib <- subset(mus, cell_type == "Dermal fibroblasts")
save_subset(mus_fib,
            file.path(out_dir, "Fibroblast", "Mus_Fib.rds"),
            "Mus Fib")

# KC (all keratinocyte + follicle + sebaceous)
mus_kc_types <- c("IFE suprabasal keratinocytes", "IFE basal keratinocytes",
                  "uHF basal keratinocytes",     "uHF suprabasal keratinocytes",
                  "Hair follicle stem cells",     "Proliferating keratinocytes",
                  "Sebaceous glands")
mus_kc <- subset(mus, cell_type %in% mus_kc_types)
save_subset(mus_kc,
            file.path(out_dir, "Keratinocyte", "Mus_KC.rds"),
            "Mus KC")

# =====================================================
# ACOMYS
# =====================================================
cat("\n=== Acomys ===\n")
aco <- readRDS(file.path(base_dir, "Acomys_data", "a.all.final.rds"))

aco_fib <- subset(aco, cell_type == "Dermal fibroblasts")
save_subset(aco_fib,
            file.path(out_dir, "Fibroblast", "Acomys_Fib.rds"),
            "Acomys Fib")

aco_kc_types <- c("IFE suprabasal keratinocytes", "IFE basal keratinocytes",
                  "uHF basal keratinocytes",     "uHF suprabasal keratinocytes",
                  "Hair follicle stem cells",     "Proliferating keratinocytes",
                  "Sebaceous glands")
aco_kc <- subset(aco, cell_type %in% aco_kc_types)
save_subset(aco_kc,
            file.path(out_dir, "Keratinocyte", "Acomys_KC.rds"),
            "Acomys KC")

cat("\n=== DONE ===\n")
cat("Output at:", out_dir, "\n")
