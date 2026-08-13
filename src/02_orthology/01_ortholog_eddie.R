# ============================================================
# 02 -- Build 4-species ortholog map -> human symbols
# Uses gprofiler2 (Ensembl Compara backend, needs internet)
# Fallback: uppercase-symbol matching if API is unreachable.
# Output: ortholog/ortholog_map.tsv
#         (columns: gene_human | gene_reindeer | gene_mouse | gene_acomys)
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  if (!requireNamespace("gprofiler2", quietly = TRUE))
    install.packages("gprofiler2", repos = "https://cloud.r-project.org")
  library(gprofiler2)
})

base <- Sys.getenv("CROSS_ROOT",
                   "/exports/eddie/scratch/YOURUSER/CrossSpecies")
sub_dir <- file.path(base, "subset")
out_dir <- file.path(base, "ortholog")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# We only need gene lists per species; use one representative subset per species
species_files <- list(
  reindeer = file.path(sub_dir, "Fibroblast", "Reindeer_Fib.rds"),
  human    = file.path(sub_dir, "Fibroblast", "Human_Fib.rds"),
  mouse    = file.path(sub_dir, "Fibroblast", "Mus_Fib.rds"),
  acomys   = file.path(sub_dir, "Fibroblast", "Acomys_Fib.rds")
)

# gprofiler2 species codes (Ensembl short names)
gpro_code <- list(
  reindeer = "btaurus",     # Reindeer aligned to bovine reference
  human    = "hsapiens",
  mouse    = "mmusculus",
  acomys   = "acahirinus"
)

gene_lists <- lapply(species_files, function(p) {
  cat("Reading", basename(p), "...\n")
  rownames(readRDS(p))
})

cat("\nGene counts by species:\n")
for (s in names(gene_lists))
  cat(sprintf("  %-10s %d\n", s, length(gene_lists[[s]])))

# --- Query 3 non-human species -> human orthologs ---
cat("\nQuerying gprofiler2 (needs internet)...\n")

query_orth <- function(genes, from_code) {
  cat("  ", from_code, "-> hsapiens\n")
  res <- tryCatch(
    gorth(query           = genes,
          source_organism = from_code,
          target_organism = "hsapiens",
          mthreshold      = 1,     # keep 1:1 only
          filter_na       = TRUE),
    error = function(e) { message("  API error: ", conditionMessage(e)); NULL }
  )
  if (is.null(res) || nrow(res) == 0) return(NULL)
  res[!duplicated(res$input), c("input", "ortholog_name")]
}

reindeer_orth <- query_orth(gene_lists$reindeer, gpro_code$reindeer)
mouse_orth    <- query_orth(gene_lists$mouse,    gpro_code$mouse)
acomys_orth   <- query_orth(gene_lists$acomys,   gpro_code$acomys)

# Human is identity mapping
human_orth <- data.frame(input = gene_lists$human,
                          ortholog_name = gene_lists$human,
                          stringsAsFactors = FALSE)

# --- Build wide ortholog table keyed by human symbol ---
merge_by_human <- function(...) {
  dfs <- list(...)
  Reduce(function(a, b)
    merge(a, b, by = "human", all = FALSE), dfs)
}

# Prepare each species df with columns human, <species>
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

# Merge on human symbol; keep only 1:1:1:1
ortho <- merge_by_human(d_human, d_reindeer, d_mouse, d_acomys)
cat("\n1:1:1:1 orthologs found:", nrow(ortho), "\n")

# Save
write.table(ortho, file.path(out_dir, "ortholog_map.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Wrote:", file.path(out_dir, "ortholog_map.tsv"), "\n")

# Sanity: print lab hero-gene coverage
hero <- c("WNT2","RSPO3","GREM1","SULF1","CHRDL1","IL13RA2","COL13A1",
          "COL1A1","DCN","LUM","PDGFRA")
cat("\nHero-gene coverage:\n")
for (g in hero) {
  ok <- g %in% ortho$human
  cat(sprintf("  %-10s %s\n", g, ifelse(ok, "OK", "missing")))
}
