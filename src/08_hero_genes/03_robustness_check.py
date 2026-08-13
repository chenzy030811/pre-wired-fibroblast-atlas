"""
09 -- Robustness check: does the hero-gene-enriched cluster persist
      under different Leiden resolutions?

For each resolution in [0.15, 0.20, 0.25, 0.30, 0.35, 0.40]:
  1. Recluster the Fib scANVI latent
  2. Find cluster with highest CHRDL1 and highest RSPO3
  3. Compute Hero-gene module score for regen vs scar in that cluster
  4. Report cluster count + module log2FC + p-value

Shows the pre-wired findings are stable across a range of clustering choices.

Output:
  results/robustness_check.csv
  preview_output/PREVIEW_v2_robustness.png
"""
import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from scipy import sparse, stats

HERO_GENES = ["WNT2", "RSPO3", "GREM1", "SULF1",
              "CHRDL1", "IL13RA2", "COL13A1"]
REGEN = ["Reindeer_Antler", "Acomys"]
SCAR  = ["Reindeer_Back", "Mus", "Human"]

base    = Path("/Users/apple/Downloads/研究生毕业论文/Cross_Species")
in_path = base / "scvi_out/Fibroblast_v2_scanvi/integrated.h5ad"
results = base / "results"; results.mkdir(exist_ok=True, parents=True)
preview = Path("/Users/apple/Downloads/preview_output")

print(">>> Loading + log-norm ...")
adata = sc.read_h5ad(in_path); adata.obs_names_make_unique()
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

hero_present = [g for g in HERO_GENES if g in adata.var_names]
X_hero = adata[:, hero_present].X
if sparse.issparse(X_hero): X_hero = X_hero.toarray()
adata.obs["Hero_score"] = X_hero.mean(axis=1)

print(">>> Rebuilding neighbors on X_scANVI ...")
sc.pp.neighbors(adata, use_rep="X_scANVI")

sp_arr = adata.obs["species"].astype(str).values
score_arr = adata.obs["Hero_score"].values

rows = []
resolutions = [0.15, 0.20, 0.25, 0.30, 0.35, 0.40]

for res in resolutions:
    print(f"\n>>> Leiden res={res} ...")
    sc.tl.leiden(adata, resolution=res, key_added=f"_res_{res}",
                 flavor="igraph", n_iterations=2, directed=False)
    cl_arr = adata.obs[f"_res_{res}"].astype(str).values
    uniq   = sorted(np.unique(cl_arr), key=int)
    k = len(uniq)

    # find TOP hero-enriched cluster that ALSO has both regen and scar cells
    # (a species-specific cluster is not useful for regen-vs-scar comparison)
    candidates = []
    for cl in uniq:
        m = cl_arr == cl
        if m.sum() < 20: continue
        n_regen = int((m & np.isin(sp_arr, REGEN)).sum())
        n_scar  = int((m & np.isin(sp_arr, SCAR)).sum())
        if n_regen < 5 or n_scar < 5: continue
        candidates.append((cl, score_arr[m].mean()))
    if not candidates:
        print(f"    (no cluster with both regen and scar ≥ 5)")
        continue

    # pick top by mean hero score
    top_cl = max(candidates, key=lambda x: x[1])[0]
    m_top = cl_arr == top_cl
    r = score_arr[m_top & np.isin(sp_arr, REGEN)]
    s = score_arr[m_top & np.isin(sp_arr, SCAR)]

    u, p = stats.mannwhitneyu(r, s, alternative="two-sided")
    lfc = np.log2((r.mean()+1e-3)/(s.mean()+1e-3))

    rows.append({
        "resolution": res, "n_clusters": k,
        "top_hero_cluster": top_cl,
        "n_cells_top": int(m_top.sum()),
        "mean_hero_score": round(float(score_arr[m_top].mean()), 3),
        "n_regen": len(r), "n_scar": len(s),
        "log2fc": round(lfc, 3),
        "p_value": p,
    })
    print(f"    k={k}, top hero cluster={top_cl}, "
          f"log2FC={lfc:+.2f}, p={p:.2e}")

df = pd.DataFrame(rows)
df.to_csv(results / "robustness_check.csv", index=False)
print(f"\n    → {results}/robustness_check.csv")

# ==== plot: log2FC + p across resolutions ====
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

ax1.plot(df["resolution"], df["log2fc"], marker="o",
         color="#D93F49", lw=2, markersize=8)
ax1.axhline(0, color="grey", lw=0.5, ls="--")
ax1.set_xlabel("Leiden resolution", fontsize=10)
ax1.set_ylabel("log2FC (regen ÷ scar) in top-hero cluster", fontsize=10)
ax1.set_title("Effect size stability", fontsize=11, fontweight="bold")
for i, r in enumerate(df["resolution"]):
    ax1.annotate(f"k={df['n_clusters'].iloc[i]}",
                 (r, df["log2fc"].iloc[i]),
                 xytext=(0, 8), textcoords="offset points",
                 fontsize=8, ha="center", color="#666666")
for s_ in ["top", "right"]: ax1.spines[s_].set_visible(False)
ax1.grid(axis="y", linestyle=":", color="#D0D0D0", linewidth=0.4)

ax2.plot(df["resolution"], -np.log10(df["p_value"]), marker="s",
         color="#5090C6", lw=2, markersize=8)
ax2.axhline(-np.log10(0.05), color="grey", lw=0.5, ls="--",
            label="p=0.05")
ax2.axhline(-np.log10(0.001), color="grey", lw=0.5, ls=":",
            label="p=0.001")
ax2.set_xlabel("Leiden resolution", fontsize=10)
ax2.set_ylabel("-log10(p)", fontsize=10)
ax2.set_title("Statistical significance stability", fontsize=11, fontweight="bold")
ax2.legend(frameon=False, fontsize=8)
for s_ in ["top", "right"]: ax2.spines[s_].set_visible(False)
ax2.grid(axis="y", linestyle=":", color="#D0D0D0", linewidth=0.4)

plt.suptitle("Robustness: top-hero cluster's regen-vs-scar effect across resolutions",
             y=1.02, fontsize=13, fontweight="bold")
plt.tight_layout()
plt.savefig(preview / "PREVIEW_v2_robustness.png", dpi=400, bbox_inches="tight")
plt.close()
print(f"    → {preview}/PREVIEW_v2_robustness.png")
print("=== DONE ===")
