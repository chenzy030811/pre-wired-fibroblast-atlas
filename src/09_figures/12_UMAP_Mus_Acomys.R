# ============================================================
# Mus + Acomys UMAP -- cloud style, matched palette to Reynolds.
# Uses pre-computed UMAP coordinates in the saved Seurat objects.
# ============================================================
suppressPackageStartupMessages({
  library(Seurat); library(tidyverse); library(patchwork); library(ggrepel)
})

# Unified 13-color deep palette (kept consistent across species)
CELL_COLORS <- c(
  "IFE basal keratinocytes"      = "#C9A227",
  "IFE suprabasal keratinocytes" = "#D17216",
  "Proliferating keratinocytes"  = "#A07C20",
  "Hair follicle stem cells"     = "#B85820",
  "uHF basal keratinocytes"      = "#8B6F2F",
  "uHF suprabasal keratinocytes" = "#9C4A1C",
  "Sebaceous glands"             = "#C66055",
  "Dermal fibroblasts"           = "#E07A5F",
  "Muscle"                       = "#1F4F94",
  "Endothelial cells"            = "#2D9686",
  "Schwann cells"                = "#85316C",
  "Immune cells"                 = "#9560AE",
  "Melanocytes"                  = "#86A53C"
)

theme_umap <- function(base = 9) {
  theme_classic(base_size = base, base_family = "Helvetica") + theme(
    plot.title = element_text(size = base + 4, face = "bold", hjust = 0),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(size = base, color = "grey40"),
    axis.line  = element_line(linewidth = 0.4),
    legend.position = "none",
    plot.margin = margin(8, 8, 8, 8))
}

# Plot helper
plot_species <- function(seu, species_name, out_path) {
  emb <- Embeddings(seu, "umap")
  df  <- data.frame(UMAP1 = emb[, 1], UMAP2 = emb[, 2],
                    cell_type = as.character(seu$cell_type))
  centers <- df %>% group_by(cell_type) %>%
    summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2),
              n = n(), .groups = "drop") %>%
    filter(n >= 30)

  n_types <- length(unique(df$cell_type))
  title <- sprintf("%s - %s cells, %d cell types",
                   species_name,
                   format(nrow(df), big.mark = ","),
                   n_types)

  p <- ggplot(df, aes(UMAP1, UMAP2, color = cell_type)) +
    geom_point(size = 0.4, alpha = 0.7, stroke = 0) +
    scale_color_manual(values = CELL_COLORS, name = NULL,
                       na.value = "#888888") +
    ggrepel::geom_text_repel(data = centers, aes(label = cell_type),
                             color = "grey15", size = 3.2,
                             bg.color = "white", bg.r = 0.18,
                             fontface = "plain", segment.color = NA,
                             max.overlaps = Inf) +
    labs(title = title, x = "UMAP 1", y = "UMAP 2") +
    theme_umap(9)

  ggsave(out_path, plot = p,
         width = 8, height = 7, dpi = 300, units = "in",
         bg = "white", device = "png", type = "cairo")
  cat("Wrote:", out_path, "\n")
  invisible(p)
}

# ---- Mus ----
cat("=== Mus ===\n")
mus <- readRDS("/Users/apple/Downloads/研究生毕业论文/Mus_data/m.all.final.rds")
cat("Mus cells:", ncol(mus), "  types:", length(unique(mus$cell_type)), "\n")
plot_species(mus, "Mus musculus",
             "/Users/apple/Downloads/preview_output/PREVIEW_Mus_UMAP.png")

# ---- Acomys ----
cat("\n=== Acomys ===\n")
aco <- readRDS("/Users/apple/Downloads/研究生毕业论文/Acomys_data/a.all.final.rds")
cat("Acomys cells:", ncol(aco), "  types:", length(unique(aco$cell_type)), "\n")
plot_species(aco, "Acomys cahirinus",
             "/Users/apple/Downloads/preview_output/PREVIEW_Acomys_UMAP.png")

cat("\nDone.\n")
