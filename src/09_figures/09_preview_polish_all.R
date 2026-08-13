# Polish-all preview — regenerates QC / HVG / Elbow / DotPlot / Integration
# with the unified Cell palette (matched to the two approved figures).
suppressPackageStartupMessages({
  library(Seurat); library(tidyverse); library(patchwork); library(ggrepel)
  library(scales)
})

base_dir <- "/Users/apple/Downloads/研究生毕业论文/reindeer_data"

# ============================================================
# Unified palette (same as approved figures)
# ============================================================
SAMPLE_COLORS <- c(Day0_Antler = "#E07A5F", Day0_Back = "#3D8DAE")
TISSUE_COLORS <- c(Velvet = "#E07A5F", Back = "#3D8DAE")

CELLTYPE_COLORS <- c(
  "Fibroblast"              = "#C66055",
  "Fibroblasts"             = "#C66055",
  "Basal Keratinocyte"      = "#A07C20",
  "Basal Keratinocytes"     = "#A07C20",
  "Suprabasal Keratinocyte" = "#AD5818",
  "Endothelium"             = "#2D9686",
  "Endothelial"             = "#2D9686",
  "Melanocyte"              = "#7C8E2C",
  "T cell"                  = "#9560AE",
  "T-cell"                  = "#9560AE",
  "CD45+IL1β+ Myeloid"      = "#5090C6",
  "Myeloid"                 = "#5090C6",
  "VSM"                     = "#1F4F94",
  "Macrophage (MΦ)"         = "#B83668",
  "Macrophages"             = "#B83668",
  "Schwann"                 = "#85316C"
)

# DotPlot gradient — warm color (matches Fibroblast / main UMAP tone)
DOTPLOT_LOW  <- "grey90"
DOTPLOT_HIGH <- "#B83668"   # deep rose; saturated, not blue-white

# Publication theme
theme_pub <- function(base = 8) {
  theme_classic(base_size = base, base_family = "Helvetica") + theme(
    plot.title    = element_text(size = base + 2, face = "bold", hjust = 0.5,
                                 margin = margin(b = 4)),
    plot.subtitle = element_text(size = base, color = "grey30", hjust = 0.5,
                                 margin = margin(b = 6)),
    axis.title    = element_text(size = base + 1),
    axis.text     = element_text(size = base, color = "black"),
    axis.line     = element_line(linewidth = 0.4, color = "black"),
    axis.ticks    = element_line(linewidth = 0.3, color = "black"),
    legend.title  = element_text(size = base, face = "bold"),
    legend.text   = element_text(size = base - 1),
    legend.key.size  = unit(0.35, "cm"),
    strip.text       = element_text(size = base + 1, face = "bold"),
    strip.background = element_blank(),
    plot.margin      = margin(6, 8, 6, 8)
  )
}
theme_umap <- function(base = 8) {
  theme_pub(base) + theme(
    axis.text = element_blank(), axis.ticks = element_blank(),
    axis.title = element_text(size = base, color = "grey40"))
}

save_png <- function(p, name, w, h) {
  ggsave(paste0("/tmp/", name, ".png"), plot = p, width = w, height = h,
         dpi = 300, units = "in", bg = "white", device = "png", type = "cairo")
}

# ============================================================
# 1. DotPlot (paper-aligned annotation) — new gradient
# ============================================================
cat("Rendering DotPlot...\n")
seu <- readRDS(file.path(base_dir, "Reindeer_Day0_annotated.rds"))
Idents(seu) <- seu$celltype_paper

key_markers <- c(
  "KRT5","KRT14","KRT15","LHX2","KRTDAP","AQP3","CST6",
  "KRT35","KRT85","DLX3","COL11A1","SFRP2","DPEP1",
  "TCF21","PCOLCE2","MFAP5","MKX","ASPN","COCH",
  "KERA","MATN4","PODN","VWF","PLVAP","SELE",
  "PROX1","CCL21","MMRN1","RGS5","MYH11","DES",
  "MLANA","TYRP1","SOX10","HMGCS2","ACSBG1","SLPI",
  "CD3E","CD8A","C1QA","C1QB","FOLR2",
  "CD1B","BOLA-DYA","CPA3","TPSB2","MS4A2"
)

p_dot <- DotPlot(seu, features = unique(key_markers),
                 cols = c(DOTPLOT_LOW, DOTPLOT_HIGH),
                 dot.scale = 4, dot.min = 0.02) +
  labs(x = NULL, y = NULL, title = "Canonical marker expression by cell type") +
  theme_pub(8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7, face = "italic"),
    axis.text.y = element_text(size = 8),
    legend.position = "right",
    panel.grid.major = element_line(color = "grey92", linewidth = 0.2)
  ) +
  guides(size  = guide_legend(title = "% expressing", order = 1),
         color = guide_colorbar(title = "Avg expr (z)", order = 2,
                                barwidth = 0.4, barheight = 4))
save_png(p_dot, "PREVIEW_DotPlot_celltype", 13, 5.5)

# ============================================================
# 2. Integration comparison — sample colors + Velvet/Back labels
# ============================================================
cat("Rendering Integration comparison...\n")
seu_n <- readRDS(file.path(base_dir, "Reindeer_Day0_normalized.rds"))

# Add tissue alias
seu_n$tissue <- ifelse(seu_n$orig.ident == "Day0_Antler", "Velvet", "Back")
seu_n$tissue <- factor(seu_n$tissue, levels = c("Velvet", "Back"))

make_int_panel <- function(red_name, title) {
  emb <- Embeddings(seu_n, red_name)[, 1:2]
  df  <- data.frame(UMAP1 = emb[, 1], UMAP2 = emb[, 2],
                    sample = seu_n$tissue)
  ggplot(df, aes(UMAP1, UMAP2, color = sample)) +
    geom_point(size = 0.3, alpha = 0.7) +
    scale_color_manual(values = TISSUE_COLORS, name = NULL) +
    labs(title = title, x = "UMAP 1", y = "UMAP 2") +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    theme_umap(8)
}

reds <- list("umap.pca"             = "No integration (PCA)",
             "umap.harmony"         = "Harmony",
             "umap.integrated.rpca" = "RPCA",
             "umap.integrated.cca"  = "CCA")
panels <- lapply(names(reds), function(r) make_int_panel(r, reds[[r]]))
fig_int <- (panels[[1]] | panels[[2]]) / (panels[[3]] | panels[[4]]) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Day 0 - integration method comparison",
    subtitle = "UMAP colored by sample (Velvet vs Back). Better integration = uniform color mixing.",
    theme    = theme(
      plot.title    = element_text(size = 12, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 8, color = "grey30", hjust = 0))) &
  theme(legend.position = "bottom")
save_png(fig_int, "PREVIEW_Integration_comparison", 9, 8)

# ============================================================
# 3. ElbowPlot — verify
# ============================================================
cat("Rendering ElbowPlot...\n")
sdev <- Stdev(seu_n, "pca")[1:30]
elbow_df <- data.frame(PC = seq_along(sdev), SD = sdev)
p_elbow <- ggplot(elbow_df, aes(PC, SD)) +
  geom_line(color = "grey60", linewidth = 0.4) +
  geom_point(color = TISSUE_COLORS[["Back"]], size = 1.8) +
  geom_vline(xintercept = 10, linetype = "dashed",
             color = TISSUE_COLORS[["Velvet"]], linewidth = 0.5) +
  annotate("text", x = 10.5, y = max(sdev) * 0.92,
           label = "Selected: 10 PCs", hjust = 0, size = 3,
           color = TISSUE_COLORS[["Velvet"]], fontface = "bold") +
  scale_x_continuous(breaks = seq(0, 30, 5)) +
  labs(x = "Principal component", y = "Standard deviation",
       title = "PCA elbow") + theme_pub(8)
save_png(p_elbow, "PREVIEW_ElbowPlot", 6, 4)

# ============================================================
# 4. HVG plot — verify
# ============================================================
cat("Rendering HVG plot...\n")
top10 <- head(VariableFeatures(seu_n), 10)
hvf_df <- HVFInfo(seu_n) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  mutate(
    variable = gene %in% VariableFeatures(seu_n),
    label    = ifelse(gene %in% top10, gene, NA)
  ) %>%
  filter(mean > 0)

p_hvg <- ggplot(hvf_df, aes(x = mean, y = variance.standardized,
                            color = variable)) +
  geom_point(size = 0.5, alpha = 0.8) +
  scale_color_manual(values = c(`TRUE`  = TISSUE_COLORS[["Velvet"]],
                                `FALSE` = "grey55"),
                     labels = c(`TRUE`  = sprintf("Variable (n = %d)",
                                                  sum(hvf_df$variable)),
                                `FALSE` = sprintf("Non-variable (n = %d)",
                                                   sum(!hvf_df$variable))),
                     breaks = c("TRUE","FALSE"), name = NULL) +
  scale_x_log10(labels = label_log()) +
  ggrepel::geom_text_repel(aes(label = label), color = "black", size = 2.6,
                           bg.color = "white", bg.r = 0.15,
                           min.segment.length = 0, max.overlaps = 20,
                           segment.color = "grey40", segment.size = 0.25,
                           na.rm = TRUE) +
  labs(x = "Average expression (log10)", y = "Standardized variance",
       title = "Highly variable genes (vst, top 2,000)") +
  theme_pub(8) +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = alpha("white", 0.8), color = NA))
suppressWarnings(save_png(p_hvg, "PREVIEW_HVG_plot", 7, 4.5))

# ============================================================
# 5. QC report — verify
# ============================================================
cat("Rendering QC report...\n")
seu_post <- readRDS(file.path(base_dir, "Reindeer_Day0_postQC.rds"))

# Reload pre-QC quickly
samples <- c("Day0_Antler", "Day0_Back")
MT_PATTERN <- "^MT-|^mt-|^Mt-|^ND[1-6]$|^ND4L$|^COX[1-3]$|^ATP[68]$|^CYTB$"
seu_list <- lapply(samples, function(s) {
  CreateSeuratObject(Read10X(file.path(base_dir, s)),
                     project = s, min.cells = 3, min.features = 200)
})
seu_pre <- merge(seu_list[[1]], seu_list[[2]],
                 add.cell.ids = samples, project = "Reindeer_Day0")
seu_pre[["percent.mt"]] <- PercentageFeatureSet(seu_pre, pattern = MT_PATTERN)

TH_FEATURE_LOW  <- 800;  TH_FEATURE_HIGH <- 4500
TH_COUNT_LOW    <- 1500; TH_COUNT_HIGH   <- 15000
TH_MT           <- 5

qc_violin <- function(df, metric, ylab, low = NULL, high = NULL, y_cap = NULL,
                      title = NULL) {
  p <- ggplot(df, aes(x = orig.ident, y = .data[[metric]], fill = orig.ident)) +
    geom_violin(scale = "width", width = 0.85, color = "grey25",
                linewidth = 0.35, alpha = 0.85) +
    geom_boxplot(width = 0.12, outlier.shape = NA, color = "grey15",
                 fill = "white", linewidth = 0.3) +
    scale_fill_manual(values = SAMPLE_COLORS, guide = "none") +
    scale_x_discrete(labels = c(Day0_Antler = "Velvet", Day0_Back = "Back")) +
    labs(x = NULL, y = ylab, title = title) +
    theme_pub(8) + theme(axis.title.x = element_blank())
  if (!is.null(low))  p <- p + geom_hline(yintercept = low, linetype = "dashed",
                                          color = "grey25", linewidth = 0.35)
  if (!is.null(high)) p <- p + geom_hline(yintercept = high, linetype = "dashed",
                                          color = "grey25", linewidth = 0.35)
  if (!is.null(y_cap)) p <- p + coord_cartesian(ylim = c(0, y_cap))
  p
}

pre_df  <- seu_pre@meta.data
post_df <- seu_post@meta.data
ycap_f <- quantile(pre_df$nFeature_RNA, 0.995)
ycap_c <- quantile(pre_df$nCount_RNA,   0.995)
ycap_m <- max(15, TH_MT * 2)

p1a <- qc_violin(pre_df,  "nFeature_RNA", "Genes / cell",
                 TH_FEATURE_LOW, TH_FEATURE_HIGH, ycap_f, "Before QC")
p1b <- qc_violin(pre_df,  "nCount_RNA",   "UMIs / cell",
                 TH_COUNT_LOW,   TH_COUNT_HIGH,   ycap_c)
p1c <- qc_violin(pre_df,  "percent.mt",   "Mitochondrial %",
                 NULL, TH_MT, ycap_m)
p2a <- qc_violin(post_df, "nFeature_RNA", "Genes / cell", title = "After QC")
p2b <- qc_violin(post_df, "nCount_RNA",   "UMIs / cell")
p2c <- qc_violin(post_df, "percent.mt",   "Mitochondrial %")

hdr <- sprintf(
  "Day 0 reindeer scRNA-seq | thresholds: nFeature %d-%d, nCount %d-%d, %%mt < %d%%
Before QC: %d cells (Velvet %d / Back %d)   After QC: %d cells (Velvet %d / Back %d, %.1f%% retained)",
  TH_FEATURE_LOW, TH_FEATURE_HIGH, TH_COUNT_LOW, TH_COUNT_HIGH, TH_MT,
  ncol(seu_pre),  sum(seu_pre$orig.ident == "Day0_Antler"),
  sum(seu_pre$orig.ident == "Day0_Back"),
  ncol(seu_post), sum(seu_post$orig.ident == "Day0_Antler"),
  sum(seu_post$orig.ident == "Day0_Back"),
  100 * ncol(seu_post) / ncol(seu_pre)
)

fig_qc <- (p1a | p1b | p1c) / (p2a | p2b | p2c) +
  plot_annotation(
    title    = "Day 0 reindeer scRNA-seq - quality control",
    subtitle = hdr,
    theme = theme(
      plot.title    = element_text(size = 11, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 7.5, color = "grey30", hjust = 0,
                                   lineheight = 1.3)))
save_png(fig_qc, "PREVIEW_QC_report_full", 10, 6.5)

# ============================================================
# 6. Integration comparison SPLIT (the one PI asked to keep)
# Replaces rainbow defaults + Day0_Antler labels + UMAP_1/2 axes
# with unified palette + Velvet/Back + UMAP 1/2.
# ============================================================
cat("Rendering Integration comparison split...\n")
seu_a <- readRDS(file.path(base_dir, "Reindeer_Day0_annotated.rds"))
seu_a$tissue <- factor(
  ifelse(seu_a$orig.ident == "Day0_Antler", "Velvet", "Back"),
  levels = c("Velvet", "Back"))

# Reuse paper cell-type palette so the split rows are visually consistent
# with the main UMAP figure. Cluster colors come from celltype_paper.
ct <- as.character(seu_a$celltype_paper)
unique_ct <- unique(ct)
SPLIT_PALETTE <- CELLTYPE_COLORS[intersect(unique_ct, names(CELLTYPE_COLORS))]
miss <- setdiff(unique_ct, names(SPLIT_PALETTE))
if (length(miss))
  SPLIT_PALETTE <- c(SPLIT_PALETTE, setNames(rep("#888888", length(miss)), miss))

make_split_panel <- function(red_name, title) {
  emb <- Embeddings(seu_a, red_name)[, 1:2]
  df  <- data.frame(UMAP1   = emb[, 1], UMAP2 = emb[, 2],
                    celltype = seu_a$celltype_paper,
                    tissue  = seu_a$tissue)
  ggplot(df, aes(UMAP1, UMAP2, color = celltype)) +
    geom_point(size = 0.55, alpha = 0.9, stroke = 0) +
    scale_color_manual(values = SPLIT_PALETTE, name = NULL) +
    facet_wrap(~ tissue, ncol = 2) +
    labs(title = title, x = "UMAP 1", y = "UMAP 2") +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    theme_umap(8) +
    theme(strip.text  = element_text(size = 9, face = "bold"),
          panel.border = element_rect(color = "grey80", fill = NA,
                                      linewidth = 0.3))
}

p_pca_s <- make_split_panel("umap.pca",             "No integration (PCA)")
p_har_s <- make_split_panel("umap.harmony",         "Harmony")
p_rpc_s <- make_split_panel("umap.integrated.rpca", "RPCA")
p_cca_s <- make_split_panel("umap.integrated.cca",  "CCA")

fig_int_split <- (p_pca_s / p_har_s / p_rpc_s / p_cca_s) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Day 0 - Integration comparison (split by sample)",
    subtitle = "Each row is one method. Cell-type structures should look similar across Velvet and Back if integration works.",
    theme = theme(
      plot.title    = element_text(size = 12, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 8, color = "grey30", hjust = 0))) &
  theme(legend.position = "bottom")

save_png(fig_int_split, "PREVIEW_Integration_comparison_split", 11, 16)

# ============================================================
# 7. UMAP by sample (post-integration) -- optimized
# Adds: cell counts in legend, interpretation subtitle, equal aspect,
# bigger points. Uses the annotated RDS (final integrated UMAP).
# ============================================================
cat("Rendering UMAP by sample (optimized)...\n")

n_velvet <- sum(seu_a$tissue == "Velvet")
n_back   <- sum(seu_a$tissue == "Back")

emb_v <- Embeddings(seu_a, "umap")
df_smp <- data.frame(
  UMAP1 = emb_v[, 1], UMAP2 = emb_v[, 2],
  tissue = factor(seu_a$tissue,
                  levels = c("Velvet", "Back"),
                  labels = c(sprintf("Velvet (n = %s)", format(n_velvet, big.mark = ",")),
                             sprintf("Back (n = %s)",   format(n_back,   big.mark = ","))))
)
# Match named palette to the new factor labels
SMP_PALETTE <- setNames(c(TISSUE_COLORS[["Velvet"]], TISSUE_COLORS[["Back"]]),
                        levels(df_smp$tissue))

# Bottom-left corner UMAP arrow helper
corner_axes <- function(df, frac = 0.10, label_size = 2.6, gap = 0.06) {
  xr <- range(df$UMAP1); yr <- range(df$UMAP2)
  xs <- xr[1] - diff(xr) * gap;  ys <- yr[1] - diff(yr) * gap
  xl <- xs + diff(xr) * frac;    yl <- ys + diff(yr) * frac
  list(
    annotate("segment", x = xs, xend = xl, y = ys, yend = ys,
             arrow = arrow(length = unit(0.08, "cm"), type = "closed"),
             linewidth = 0.4),
    annotate("segment", x = xs, xend = xs, y = ys, yend = yl,
             arrow = arrow(length = unit(0.08, "cm"), type = "closed"),
             linewidth = 0.4),
    annotate("text", x = xs + diff(xr) * frac * 0.5,
             y = ys - diff(yr) * 0.04, label = "UMAP 1",
             size = label_size, hjust = 0.5, vjust = 1),
    annotate("text", x = xs - diff(xr) * 0.03,
             y = ys + diff(yr) * frac * 0.5, label = "UMAP 2",
             size = label_size, hjust = 0.5, vjust = 0.5, angle = 90)
  )
}
expand_low <- function(rng, gap = 0.18) {
  d <- diff(rng); c(rng[1] - d * gap, rng[2])
}

p_smp_opt <- ggplot(df_smp, aes(UMAP1, UMAP2, color = tissue)) +
  geom_point(size = 0.5, alpha = 0.8, stroke = 0) +
  scale_color_manual(values = SMP_PALETTE, name = NULL) +
  corner_axes(df_smp) +
  coord_fixed(xlim = expand_low(range(df_smp$UMAP1)),
              ylim = expand_low(range(df_smp$UMAP2)),
              clip = "off") +
  labs(title    = "Sample distribution (post-integration)",
       subtitle = "Velvet and Back cells co-mingle within each cluster - Harmony integration successful") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme_classic(base_size = 9, base_family = "Helvetica") +
  theme(
    plot.title       = element_text(size = 12, face = "bold", hjust = 0),
    plot.subtitle    = element_text(size = 9, color = "grey30", hjust = 0,
                                    margin = margin(b = 8)),
    axis.text  = element_blank(), axis.ticks = element_blank(),
    axis.title = element_blank(), axis.line  = element_blank(),
    legend.position  = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = alpha("white", 0.85), color = NA),
    legend.text      = element_text(size = 9),
    legend.key.size  = unit(0.4, "cm"),
    plot.margin      = margin(8, 12, 8, 12))

save_png(p_smp_opt, "PREVIEW_UMAP_by_sample", 7, 6)

cat("\nAll 7 previews written to /tmp/PREVIEW_*.png\n")
