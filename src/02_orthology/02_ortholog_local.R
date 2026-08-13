# ============================================================
# 02 -- Build 4-species ortholog map -> human symbols  (LOCAL Mac version)
# Reindeer / Mouse : queried via gprofiler2 (Ensembl Compara)
# Acomys           : direct uppercase-symbol match (gprofiler doesn't have it)
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  if (!requireNamespace("gprofiler2", quietly = TRUE))
    install.packages("gprofiler2", repos = "https://cloud.r-project.org")
  library(gprofiler2)
})

base    <- "/Users/apple/Downloads/研究生毕业论文/Cross_Species"
sub_dir <- file.path(base, "subset")
out_dir <- file.path(base, "ortholog")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

species_files <- list(
  reindeer = file.path(sub_dir, "Fibroblast", "Reindeer_Fib.rds"),
  human    = file.path(sub_dir, "Fibroblast", "Human_Fib.rds"),
  mouse    = file.path(sub_dir, "Fibroblast", "Mus_Fib.rds"),
  acomys   = file.path(sub_dir, "Fibroblast", "Acomys_Fib.rds")
)

gene_lists <- lapply(species_files, function(p) {
  cat("Reading", basename(p), "...\n")
  rownames(readRDS(p))
})

cat("\nGene counts by species:\n")
for (s in names(gene_lists))
  cat(sprintf("  %-10s %d\n", s, length(gene_lists[[s]])))

# --- Ortholog querying ---

query_orth <- function(genes, from_code) {
  cat("  gprofiler2:  ", from_code, "-> hsapiens\n")
  res <- tryCatch(
    gorth(query           = genes,
          source_organism = from_code,
          target_organism = "hsapiens",
          mthreshold      = 1,
          filter_na       = TRUE),
    error = function(e) { message("  API error: ", conditionMessage(e)); NULL }
  )
  if (is.null(res) || nrow(res) == 0) return(NULL)
  res[!duplicated(res$input), c("input", "ortholog_name")]
}

cat("\n>>> Reindeer (Ensembl bovine) -> Human\n")
reindeer_orth <- query_orth(gene_lists$reindeer, "btaurus")

cat("\n>>> Mouse -> Human\n")
mouse_orth    <- query_orth(gene_lists$mouse, "mmusculus")

cat("\n>>> Acomys -> Human (direct uppercase-symbol match, Acomys not in Ensembl)\n")
# Acomys gene names are already uppercase (e.g. COL1A1, ACCS). We accept as
# ortholog any Acomys gene whose exact uppercase symbol also exists in human.
common <- intersect(gene_lists$acomys, gene_lists$human)
acomys_orth <- data.frame(input = common,
                          ortholog_name = common,
                          stringsAsFactors = FALSE)
cat("  Direct matches:", nrow(acomys_orth), "genes\n")

# Human is identity
human_orth <- data.frame(input = gene_lists$human,
                          ortholog_name = gene_lists$human,
                          stringsAsFactors = FALSE)

# --- Merge to a 1:1:1:1 table on human symbol ---
prep <- function(df, species_col) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  d <- df[!duplicated(df$ortholog_name), ]
  data.frame(human = d$ortholog_name,
             setNames(list(d$input), species_col),
             stringsAsFactors = FALSE)
}

d_human    <- prep(human_orth,    "human")
d_reindeer <- prep(reindeer_orth, "reindeer")
d_mouse    <- prep(mouse_orth,    "mouse")
d_acomys   <- prep(acomys_orth,   "acomys")

# Only merge species that returned results
all_dfs <- Filter(Negate(is.null), list(human=d_human, reindeer=d_reindeer,
                                        mouse=d_mouse, acomys=d_acomys))
cat("\nSpecies with successful ortholog mapping:", paste(names(all_dfs), collapse=", "), "\n")

ortho <- Reduce(function(a, b) merge(a, b, by = "human", all = FALSE), all_dfs)
cat("\n1:1:1:1 orthologs across all 4 species:", nrow(ortho), "\n")

write.table(ortho, file.path(out_dir, "ortholog_map.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Wrote:", file.path(out_dir, "ortholog_map.tsv"), "\n")

# --- Sanity: lab hero-gene coverage ---
hero <- c("WNT2","RSPO3","GREM1","SULF1","CHRDL1","IL13RA2","COL13A1",
          "COL1A1","DCN","LUM","PDGFRA")
cat("\nHero-gene coverage (should mostly be OK):\n")
for (g in hero) {
  ok <- g %in% ortho$human
  cat(sprintf("  %-10s %s\n", g, ifelse(ok, "OK", "MISSING")))
}
