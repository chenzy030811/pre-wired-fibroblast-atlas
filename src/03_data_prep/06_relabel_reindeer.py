"""
03c — Split Reindeer into 2 groups based on tissue (Antler vs Back).

Input:
  h5ad/merged_Fibroblast.h5ad     (species has 4 categories)
  h5ad/merged_Keratinocyte.h5ad

Output:
  h5ad/merged_Fibroblast_v2.h5ad  (species has 5 categories)
  h5ad/merged_Keratinocyte_v2.h5ad

The 5 categories:
  Reindeer_Antler  (regen)  — velvet skin
  Reindeer_Back    (scar)   — dorsal body skin
  Acomys           (regen)
  Mus              (scar)
  Human            (scar)

Also adds a 'regen_status' column with values 'Regen' or 'Scar'.
"""
import anndata as ad
import pandas as pd
from pathlib import Path

base = Path("/Users/apple/Downloads/研究生毕业论文/Cross_Species/h5ad")

for comp in ["Fibroblast", "Keratinocyte"]:
    in_path  = base / f"merged_{comp}.h5ad"
    out_path = base / f"merged_{comp}_v2.h5ad"
    print(f"\n=== {comp} ===")
    adata = ad.read_h5ad(in_path)
    adata.obs_names_make_unique()
    print(f"    Loaded {adata.n_obs} cells")

    # ---- rename species based on orig.ident for Reindeer ----
    old_species = adata.obs["species"].astype(str).values.copy()
    new_species = old_species.copy()

    reindeer_mask = old_species == "Reindeer"
    orig = adata.obs["orig.ident"].astype(str).values
    n_reindeer = reindeer_mask.sum()

    # map orig.ident → new label
    antler_mask = reindeer_mask & (pd.Series(orig).str.contains("Antler", case=False).values)
    back_mask   = reindeer_mask & (pd.Series(orig).str.contains("Back", case=False).values)

    new_species[antler_mask] = "Reindeer_Antler"
    new_species[back_mask]   = "Reindeer_Back"

    # sanity: every reindeer cell should be labelled
    n_labelled = antler_mask.sum() + back_mask.sum()
    assert n_labelled == n_reindeer, (
        f"Missed some Reindeer cells: {n_labelled}/{n_reindeer} labelled. "
        f"Unique orig.ident values for Reindeer: "
        f"{pd.Series(orig[reindeer_mask]).unique()}"
    )

    adata.obs["species"] = pd.Categorical(new_species)
    print(f"    Reindeer split: {antler_mask.sum()} Antler + {back_mask.sum()} Back")
    print(f"    New species distribution:")
    print(adata.obs["species"].value_counts())

    # ---- add regen_status column ----
    REGEN = {"Reindeer_Antler", "Acomys"}
    SCAR  = {"Reindeer_Back", "Mus", "Human"}
    status = ["Regen" if s in REGEN else "Scar" for s in new_species]
    adata.obs["regen_status"] = pd.Categorical(status, categories=["Regen", "Scar"])
    print(f"\n    Regen vs Scar:")
    print(adata.obs["regen_status"].value_counts())

    # save
    adata.write_h5ad(out_path)
    print(f"    → wrote {out_path}")

print("\n=== DONE ===")
