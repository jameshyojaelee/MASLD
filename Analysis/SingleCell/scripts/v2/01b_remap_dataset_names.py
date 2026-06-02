"""Patch hepatocyte_atlas_v2_annotated.h5ad: remap dataset '0'-'6' -> real names.

01_integrate_atlas_v2 used `ad.concat(label='dataset')` which overwrote the
real dataset names with positional integers. scVI/scANVI training used the
integer-coded categorical as batch_key — that's fine, but downstream joins
on (sample, dataset) against v1 donor_metadata.tsv (which uses real names)
silently miss-match to 0 rows. Remap here without retraining.

Mapping (alphabetical sort of per_dataset/*_v2.h5ad globbed paths):
  0 -> GSE136103
  1 -> GSE174748
  2 -> GSE185477
  3 -> GSE189600
  4 -> GSE202379
  5 -> GSE244832
  6 -> Liver_Atlas
"""
from pathlib import Path
import anndata as ad

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
ATLAS = ROOT / "Analysis/SingleCell/results_gpu_v2_phase05/atlas/hepatocyte_atlas_v2_annotated.h5ad"

MAPPING = {
    "0": "GSE136103",
    "1": "GSE174748",
    "2": "GSE185477",
    "3": "GSE189600",
    "4": "GSE202379",
    "5": "GSE244832",
    "6": "Liver_Atlas",
}

print(f"[load] {ATLAS}")
a = ad.read_h5ad(ATLAS)
print(f"  shape: {a.shape}")
print(f"  current dataset categories: {a.obs['dataset'].cat.categories.tolist()}")
print(f"  current value_counts:\n{a.obs['dataset'].value_counts()}")

new_cats = [MAPPING[c] for c in a.obs["dataset"].cat.categories]
a.obs["dataset"] = a.obs["dataset"].cat.rename_categories(new_cats)
print(f"  new dataset categories: {a.obs['dataset'].cat.categories.tolist()}")
print(f"  new value_counts:\n{a.obs['dataset'].value_counts()}")

print(f"[write] {ATLAS}")
a.write_h5ad(ATLAS, compression="gzip")
print("DONE")
