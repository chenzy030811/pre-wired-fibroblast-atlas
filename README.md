# Pre-wired Regenerative Fibroblast Atlas

Cross-species single-cell RNA-seq atlas testing whether the seven-gene
"pre-wired" regenerative fibroblast signature first observed in *Acomys
cahirinus* (Cairo spiny mouse) is conserved in other adult mammals that
heal skin without scarring.

MSc Bioinformatics dissertation code — University of Edinburgh, 2025–26.
Ferreira laboratory.

---

## ⚠️  Running on your own machine — path setup

Several R scripts in `src/01_qc/`, `src/02_orthology/` and `src/03_data_prep/`
contain hard-coded absolute paths pointing to the author's local dataset
directory (e.g. `/Users/apple/Downloads/研究生毕业论文/human_data`). Before
running any script on a different machine, replace these paths with the
location of the relevant raw dataset on your own filesystem. A single
find-and-replace across the `src/` tree is enough:

```bash
# example — from the repo root
grep -rl "/Users/apple/Downloads" src/ | \
  xargs sed -i '' 's|/Users/apple/Downloads/研究生毕业论文|/your/local/path|g'
```

The Python scripts under `src/03_data_prep/` onwards use relative paths
resolved from the repo root and do not need editing. Path-editing is a
known limitation of this dissertation-scope release; a follow-up refactor
will move all paths into a single YAML config file.

---


---

## Project scope

Five biological groups, ten datasets, one integrated atlas:

| Group | Species | Site | Phenotype |
|-------|---------|------|-----------|
| Mus | *Mus musculus* | dorsal skin | scarrer |
| Acomys | *Acomys cahirinus* | dorsal skin | regenerator |
| Human | *Homo sapiens* | healthy skin | scarrer |
| Reindeer velvet | *Rangifer tarandus* | velvet antler | regenerator |
| Reindeer back | *Rangifer tarandus* | dorsal skin | scarrer |

Seven hero genes tested a priori: **WNT2, RSPO3, GREM1, SULF1, CHRDL1,
IL13RA2, COL13A1**.

---

## Repository layout

```
src/
├── 01_qc/            Per-dataset quality control (R and Python)
├── 02_orthology/     Cross-species 1:1:1:1 ortholog mapping (gprofiler2)
├── 03_data_prep/     Subsetting, h5ad building, reindeer relabelling
├── 04_integration/   scVI and scANVI training
├── 05_benchmark/     scIB integration benchmark (compartment-specific)
├── 06_clustering/    Leiden reclustering and inspection
├── 07_markers/       Marker-gene analysis (Wilcoxon)
├── 08_hero_genes/    Module score, subtype-matched testing, robustness
└── 09_figures/       Plotting scripts for every figure in the thesis
```

Each stage directory contains its own `README.md` describing the scripts
inside.

---

## Reproducibility

- All Python code runs on Python 3.11 with the pinned versions in
  `environment.yml`; recreate with:
  ```bash
  conda env create -f environment.yml
  conda activate prewired
  ```
- All R code runs on R 4.4 with Seurat v5, gprofiler2, and standard
  Bioconductor dependencies (see `environment.yml`).
- Random seeds are fixed at 0 for every stochastic step (HVG selection,
  scVI/scANVI training, Leiden, UMAP, Scrublet, permutation null).
- GPU-dependent scripts (04_integration/) were run on the University of
  Edinburgh Eddie3 HPC cluster (SGE, GPU queue, CUDA 12.1). A single
  NVIDIA A100 completes the fibroblast scANVI training in ≈ 25 min.

---

## How to reproduce a figure from scratch

1. Download raw data (accessions in `docs/data_sources.md`).
2. Run `src/01_qc/` for each dataset.
3. Build the merged h5ad: `src/03_data_prep/05_build_h5ad.py`.
4. Compute the ortholog table: `src/02_orthology/02_ortholog_local.R`.
5. Train the integrators: `src/04_integration/` (GPU recommended).
6. Score the benchmark: `src/05_benchmark/02_scib_benchmark_v2.py`.
7. Recluster: `src/06_clustering/01_reclustering_fib_v2.py`.
8. Hero-gene test: `src/08_hero_genes/02_subtype_matched_hero.py`.
9. Plot: `src/09_figures/`.

The end-to-end wall-clock on a Slurm/SGE node with one A100 GPU is
≈ 3 hours from raw h5ad to final figures.

---

## Data availability

Raw sequencing accessions are listed in `docs/data_sources.md`. The
1:1:1:1 ortholog table used throughout is provided at
`data/orthologs_1to1to1to1.csv` (10 554 rows).

Intermediate AnnData objects (post-QC, post-integration, post-
clustering) are hosted at Zenodo — DOI to be added on publication of the
accompanying manuscript.

---

## Citation

If you use this code or the ortholog table, please cite the
dissertation (details on the cover sheet) and the underlying
Ferreira-lab Acomys atlas (manuscript in preparation).

---

## Licence

MIT (see `LICENSE`).

---

## Contact

Student email — see the dissertation cover sheet.
