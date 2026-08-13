#!/usr/bin/env bash
# End-to-end driver — runs the full pipeline in canonical order.
# Assumes raw datasets are already downloaded (see docs/data_sources.md).
set -euo pipefail
export PYTHONHASHSEED=0

echo "=== 01 · Per-dataset QC ===";           Rscript src/01_qc/reindeer_QC.R
echo "=== 02 · Ortholog mapping ===";         Rscript src/02_orthology/02_ortholog_local.R
echo "=== 03 · Merge into single h5ad ===";   python src/03_data_prep/05_build_h5ad.py
echo "=== 04 · Integration (GPU) ===";        python src/04_integration/03_train_scvi_v2.py
                                              python src/04_integration/04_train_scanvi_v2.py
echo "=== 05 · scIB benchmark ===";           python src/05_benchmark/02_scib_benchmark_v2.py
echo "=== 06 · Leiden reclustering ===";      python src/06_clustering/01_reclustering_fib_v2.py
                                              python src/06_clustering/02_reclustering_kc_scanvi.py
echo "=== 07 · Markers ===";                  python src/07_markers/03_fib_v2_markers.py
                                              python src/07_markers/04_kc_v2_markers.py
echo "=== 08 · Hero-gene testing ===";        python src/08_hero_genes/02_subtype_matched_hero.py
                                              python src/08_hero_genes/03_robustness_check.py
echo "=== 09 · Figures ===";                  python src/09_figures/07_plot_v2_umaps_labeled.py
                                              python src/09_figures/01_compose_main_figures.py

echo "All done — figures written to figures/"
