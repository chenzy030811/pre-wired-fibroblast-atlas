# ============================================================
# Reynolds 2021 UMAP -- cloud-style, paper-quality
# Pipeline: Normalize -> HVG -> PCA -> Harmony (by donor) -> UMAP
# Colored by broad cell type (40 fine -> 13 broad categories)
# ============================================================
suppressPackageStartupMessages({
  library(Seurat); library(tidyverse); library(patchwork); library(ggrepel)
  library(harmony)
})
options(future.globals.maxSize = 16 * 1024^3)

base_dir <- "/Users/apple/Downloads/研究生毕业论文/human_data"

# ---- Load post-QC Reynolds (60K cells, 5% mt for cross-species fairness) ----
seu <- readRDS(file.path(base_dir, "Reynolds_postQC.rds"))
cat("Loaded:", ncol(seu), "cells from", length(unique(seu$sample)), "donors\n")

# ---- Map fine cell types into broad categories ----
broad_map <- function(x) {
  case_when(
    grepl("Keratinocyte", x)            ~ "Keratinocyte",
    grepl("^Fibroblast",  x)            ~ "Fibroblast",
    x %in% c("Th","Tc","Treg")          ~ "T cell",
    grepl("^NK", x)                     ~ "NK cell",
    x %in% c("Macrophage 1","Macrophage 2","Monocyte") ~ "Macrophage/Monocyte",
    grepl("DC", x) | grepl("^LC", x) | grepl("KLF10 LC", x) ~ "Dendritic / Langerhans",
    grepl("Vascular endothelium", x)    ~ "Vascular endothelium",
    grepl("Lymphatic endothelium", x)   ~ "Lymphatic endothelium",
    x == "Pericyte"                     ~ "Pericyte",
    grepl("Schwann", x)                 ~ "Schwann",
    x == "Melanocyte"                   ~ "Melanocyte",
    x == "Mast cell"                    ~ "Mast cell",
    x == "Plasma cell"                  ~ "Plasma cell",
    x == "ILC"                          ~ "T cell",       # innate lymphoid
    TRUE                                ~ "Other"
  )
}
seu$broad_type <- broad_map(seu$cell_type)
cat("\nBroad cell types:\n"); print(sort(table(seu$broad_type), decreasing = TRUE))

# ---- Pipeline ----
DefaultAssay(seu) <- "RNA"
seu[["RNA"]] <- as(seu[["RNA"]], "Assay5")
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$sample)

cat("\nRunning pipeline (this takes ~5-10 min)...\n")
seu <- seu %>%
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(nfeatures = 3000, verbose = FALSE) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE)

# Harmony integration by donor (sample = s1/s2/s3)
seu <- IntegrateLayers(seu, method = HarmonyIntegration,
                       orig.reduction = "pca", new.reduction = "harmony",
                       verbose = FALSE)

seu <- seu %>%
  FindNeighbors(reduction = "harmony", dims = 1:20, verbose = FALSE) %>%
  RunUMAP(reduction = "harmony", dims = 1:20,
          min.dist = 0.3, n.neighbors = 30, spread = 1,
          reduction.name = "umap", verbose = FALSE)
seu <- JoinLayers(seu)

saveRDS(seu, file.path(base_dir, "Reynolds_umap.rds"))
cat("Saved: Reynolds_umap.rds\n\n")

# ---- Plot (cloud style) ----
# 13 deep pastel colors for broad cell types
BROAD_COLORS <- c(
  "Keratinocyte"             = "#C9A227",  # deep gold
  "Fibroblast"               = "#C66055",  # deep salmon
  "T cell"                   = "#9560AE",  # deep purple
  "NK cell"                  = "#7C8E2C",  # olive
  "Macrophage/Monocyte"      = "#B83668",  # deep rose
  "Dendritic / Langerhans"   = "#5090C6",  # sky blue
  "Vascular endothelium"     = "#2D9686",  # teal
  "Lymphatic endothelium"    = "#4DAFA0",  # lighter teal
  "Pericyte"                 = "#1F4F94",  # deep blue
  "Schwann"                  = "#85316C",  # magenta
  "Melanocyte"               = "#86A53C",  # bright olive
  "Mast cell"                = "#E07A5F",  # coral
  "Plasma cell"              = "#A04035",  # brick
  "Other"                    = "#888888"
)

emb <- Embeddings(seu, "umap")
df  <- data.frame(UMAP1 = emb[, 1], UMAP2 = emb[, 2],
                  cell_type = seu$broad_type)
centers <- df %>% group_by(cell_type) %>%
  summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2),
            n = n(), .groups = "drop") %>%
  filter(n >= 50)  # don't label tiny groups

p <- ggplot(df, aes(UMAP1, UMAP2, color = cell_type)) +
  geom_point(size = 0.4, alpha = 0.7, stroke = 0) +
  scale_color_manual(values = BROAD_COLORS, name = NULL) +
  ggrepel::geom_text_repel(data = centers, aes(label = cell_type),
                           color = "grey15", size = 3.2,
                           bg.color = "white", bg.r = 0.18,
                           fontface = "plain", segment.color = NA,
                           max.overlaps = Inf) +
  labs(title = sprintf("Reynolds 2021 - %s cells, %d cell types",
                       format(ncol(seu), big.mark=","), 13),
       x = "UMAP 1", y = "UMAP 2") +
  guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  theme_classic(base_size = 9, base_family = "Helvetica") +
  theme(
    plot.title       = element_text(size = 13, face = "bold", hjust = 0),
    axis.text        = element_blank(),
    axis.ticks       = element_blank(),
    axis.title       = element_text(size = 9, color = "grey40"),
    axis.line        = element_line(linewidth = 0.4),
    legend.position  = "none",
    plot.margin      = margin(8, 8, 8, 8)
  )

out_png <- "/Users/apple/Downloads/preview_output/PREVIEW_Reynolds_UMAP.png"
ggsave(out_png, plot = p,
       width = 8, height = 7, dpi = 300, units = "in",
       bg = "white", device = "png", type = "cairo")
cat("Wrote:", out_png, "\n")
