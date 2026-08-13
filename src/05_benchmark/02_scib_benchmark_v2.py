"""
08b -- scIB benchmark v2: score scVI vs scANVI on 5-group integrated data
       for BOTH Fibroblast and Keratinocyte.

Output:
  results/scib_benchmark_v2_Fib.csv
  results/scib_benchmark_v2_KC.csv
  preview_output/PREVIEW_scib_v2_Fib.png
  preview_output/PREVIEW_scib_v2_KC.png
  Terminal: winner for each compartment.
"""
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection

base    = Path("/Users/apple/Downloads/研究生毕业论文/Cross_Species")
results = base / "results"; results.mkdir(exist_ok=True, parents=True)
preview = Path("/Users/apple/Downloads/preview_output")

for comp in ["Fibroblast", "Keratinocyte"]:
    print(f"\n{'='*60}\n{comp}  ·  scIB v2 benchmark\n{'='*60}")

    scvi_h5   = base / f"scvi_out/{comp}_v2_scvi/integrated.h5ad"
    scanvi_h5 = base / f"scvi_out/{comp}_v2_scanvi/integrated.h5ad"

    print(">>> Loading ...")
    a_sv = ad.read_h5ad(scvi_h5);   a_sv.obs_names_make_unique()
    a_sa = ad.read_h5ad(scanvi_h5); a_sa.obs_names_make_unique()

    assert a_sv.n_obs == a_sa.n_obs, f"cell count mismatch {a_sv.n_obs} vs {a_sa.n_obs}"
    print(f"    {a_sv.n_obs} cells (matched)")

    # combined adata using scVI as base, add scANVI embedding
    adata = a_sv.copy()
    adata.obsm["X_scVI"]   = a_sv.obsm["X_scVI"]
    adata.obsm["X_scANVI"] = a_sa.obsm["X_scANVI"]

    # labels: use scVI's leiden as consensus proxy
    adata.obs["label"] = adata.obs["leiden_scvi"].astype("category")
    adata.obs["batch"] = adata.obs["species"].astype("category")

    print(">>> Running scIB benchmark ...")
    bm = Benchmarker(
        adata,
        batch_key="batch",
        label_key="label",
        embedding_obsm_keys=["X_scVI", "X_scANVI"],
        n_jobs=-1,
        bio_conservation_metrics=BioConservation(
            isolated_labels=True,
            nmi_ari_cluster_labels_leiden=True,
            silhouette_label=True,
            clisi_knn=True,
        ),
        batch_correction_metrics=BatchCorrection(
            graph_connectivity=True,
            kbet_per_label=True,
            ilisi_knn=True,
            pcr_comparison=True,
        ),
    )
    bm.benchmark()

    df = bm.get_results(min_max_scale=False)
    csv_path = results / f"scib_benchmark_v2_{comp[:3]}.csv"
    df.to_csv(csv_path)
    print(f"    → {csv_path}")

    # ---- bar chart ----
    plot_df = df.T
    summary_cols = ["Bio conservation", "Batch correction", "Total"]
    metric_rows = [i for i in plot_df.index
                   if i != "Metric Type" and i not in summary_cols]
    plot_df = plot_df.loc[metric_rows]

    fig, ax = plt.subplots(figsize=(10, 5))
    xpos = np.arange(len(plot_df))
    w = 0.36
    scvi_vals   = plot_df["X_scVI"].astype(float).values
    scanvi_vals = plot_df["X_scANVI"].astype(float).values
    ax.bar(xpos - w/2, scvi_vals,   w, label="scVI",   color="#8AB2D0",
           edgecolor="black", linewidth=0.4)
    ax.bar(xpos + w/2, scanvi_vals, w, label="scANVI", color="#E28A7A",
           edgecolor="black", linewidth=0.4)
    for i, v in enumerate(scvi_vals):
        ax.text(i-w/2, v+0.015, f"{v:.2f}", ha="center", va="bottom", fontsize=8)
    for i, v in enumerate(scanvi_vals):
        ax.text(i+w/2, v+0.015, f"{v:.2f}", ha="center", va="bottom", fontsize=8)
    ax.set_xticks(xpos)
    ax.set_xticklabels(plot_df.index, rotation=25, ha="right", fontsize=9)
    ax.set_ylabel("Score (higher = better)")
    ax.set_ylim(0, 1.08)
    ax.legend(frameon=False, fontsize=10)
    for s in ["top", "right"]: ax.spines[s].set_visible(False)
    ax.grid(axis="y", linestyle=":", linewidth=0.4, color="#D0D0D0")
    ax.set_axisbelow(True)
    plt.tight_layout()
    png_path = preview / f"PREVIEW_scib_v2_{comp[:3]}.png"
    plt.savefig(png_path, dpi=400, bbox_inches="tight")
    plt.close()
    print(f"    → {png_path}")

    # ---- winner ----
    if "Total" in df.index:
        t_sv = float(df.loc["Total", "X_scVI"])
        t_sa = float(df.loc["Total", "X_scANVI"])
        winner = "scANVI" if t_sa > t_sv else "scVI"
        print(f"\n    🏆 {comp} WINNER: {winner}  (scVI={t_sv:.3f}  scANVI={t_sa:.3f})")

print("\n=== DONE ===")
