#!/usr/bin/env python3
"""Create donor/cell-type manifests for adult-liver ChromBPNet preparation."""
import csv, json, os, re
from collections import Counter, defaultdict
import anndata as ad

ROOT = os.environ.get("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
ATAC = os.path.join(ROOT, "Analysis/ATAC/Human_Multiome")
OUT = os.path.join(ROOT, "GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet")
os.makedirs(os.path.join(OUT, "barcodes"), exist_ok=True)
h5 = os.path.join(ATAC, "results/label_transfer/snapatac2_label_transferred.h5ad")
a = ad.read_h5ad(h5, backed="r")
obs = a.obs[["donor_id", "condition", "cell_type"]].copy()

aliases = {"Hepatocyte": "hepatocyte", "Stellate_Cell": "stellate",
           "Macrophage": "macrophage", "Endothelial_Cell": "endothelial",
           "Cholangiocyte": "cholangiocyte"}
obs["model_celltype"] = [aliases.get(str(x), re.sub(r"[^a-z0-9]+", "_", str(x).lower()).strip("_")) for x in obs.cell_type]
counts = obs.groupby(["model_celltype", "donor_id", "condition"], observed=True).size().reset_index(name="n_cells")
counts.to_csv(os.path.join(OUT, "cell_counts_by_donor.tsv"), sep="\t", index=False)

for (ct, donor), idx in obs.groupby(["model_celltype", "donor_id"], observed=True).groups.items():
    with open(os.path.join(OUT, "barcodes", f"{donor}.{ct}.txt"), "w") as f:
        f.write("\n".join(map(str, idx)) + "\n")

peak_dir = os.path.join(ATAC, "results/label_transfer/cell_type_peak_sets_v2")
peak_map = {"hepatocyte": "Hepatocytes_peaks.bed", "stellate": "Fibroblasts_peaks.bed",
            "macrophage": "Macrophages_peaks.bed", "endothelial": "Endothelial_cells_peaks.bed",
            "cholangiocyte": "Cholangiocytes_peaks.bed"}
manifest = []
for ct in peak_map:
    z = counts[counts.model_celltype == ct]
    n_cells, n_donors = int(z.n_cells.sum()), int(z.donor_id.nunique())
    cond = ";".join(f"{k}:{int(v)}" for k,v in z.groupby("condition").n_cells.sum().items())
    gate = n_cells >= 10000 and n_donors >= 6 and os.path.exists(os.path.join(peak_dir, peak_map[ct]))
    manifest.append([ct, n_cells, n_donors, cond, os.path.join(peak_dir, peak_map[ct]), gate,
                     "donor-held-out chromosome-fold training; pooled model first"])
with open(os.path.join(OUT, "model_manifest.tsv"), "w", newline="") as f:
    w=csv.writer(f,delimiter="\t"); w.writerow(["cell_type","n_cells","n_donors","condition_cells","peaks","preflight_pass","design"]); w.writerows(manifest)
with open(os.path.join(OUT, "contract.json"), "w") as f:
    json.dump({"status":"inputs_prepared_training_not_started", "apply_only_firewall":True,
      "fold_rule":"donor-disjoint pseudobulks and held-out chromosomes; never cell-random split",
      "primary_models":["hepatocyte","stellate"], "secondary_models":["macrophage","endothelial","cholangiocyte"],
      "blocker":"ChromBPNet training package is not installed; approval is required before creating the isolated training environment."}, f, indent=2)
print(f"wrote {OUT}")
