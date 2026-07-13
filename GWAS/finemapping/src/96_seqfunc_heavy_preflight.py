#!/usr/bin/env python3
"""Preflight the compute-heavy seqfunc extensions without mutating canonical results."""
import csv, glob, gzip, hashlib, json, os, shutil
from collections import Counter, defaultdict

ROOT = os.environ.get("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUT = os.path.join(ROOT, "GWAS/finemapping/results/seqfunc/heavy_preflight")
os.makedirs(OUT, exist_ok=True)

def present(path):
    return {"path": os.path.relpath(path, ROOT), "present": os.path.exists(path),
            "bytes": os.path.getsize(path) if os.path.isfile(path) else None}

checks = []
atac = os.path.join(ROOT, "Analysis/ATAC/Human_Multiome")
checks += [present(os.path.join(atac, "results/label_transfer/snapatac2_label_transferred.h5ad")),
           present(os.path.join(atac, "metadata/donor_metadata_curated.tsv")),
           present(os.path.join(ROOT, "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv")),
           present(os.path.join(ROOT, "GWAS/finemapping/results/combined_finemapping.csv")),
           present(os.path.join(ROOT, "GWAS/finemapping/results/seqfunc/variant_substrate_hg38.tsv"))]

fragments = sorted(glob.glob(os.path.join(atac, "results/fragments/*_fragments.tsv.gz")))
peaks = sorted(glob.glob(os.path.join(atac, "results/label_transfer/cell_type_peak_sets_v2/*_peaks.bed")))
sj = sorted(glob.glob(os.path.join(ROOT, "RNA-seq/Human/Patient_Cohorts/results/*/alignments/star/*/*SJ.out.tab")))
ld = sorted(glob.glob(os.path.join(ROOT, "GWAS/finemapping/data/ld_ref/topld_*/*/*/*.ld")))

rows = [
    ("adult_liver_fragments", len(fragments), len(fragments) == 18, "18 donor fragment files expected"),
    ("adult_liver_celltype_peaks", len(peaks), len(peaks) >= 5, "primary liver lineages expected"),
    ("star_junction_files", len(sj), len(sj) >= 800, "cohort-scale junction substrate"),
    ("local_ancestry_ld_matrices", len(ld), len(ld) > 0, "LD matrices permit LD-aware pairs, not observed haplotypes"),
]
with open(os.path.join(OUT, "resource_inventory.tsv"), "w", newline="") as f:
    w = csv.writer(f, delimiter="\t"); w.writerow(["resource", "n", "gate_pass", "interpretation"]); w.writerows(rows)

contract = {
  "status": "preflight_complete", "apply_only_firewall": True,
  "canonical_outputs_mutated": False,
  "input_checks": checks,
  "resources": {r[0]: {"n": r[1], "gate_pass": r[2]} for r in rows},
  "environment_gates": {
    "chrombpnet_training_package": "MISSING; chrombpnet_vs contains the scorer/TensorFlow but not the training package",
    "leafcutter": "MISSING; junction inventory/pilot can proceed, production clustering requires an isolated environment",
    "phased_reference_genotypes": "MISSING locally; .ld/.bim windows are not sufficient to reconstruct observed haplotypes"
  },
  "hard_stops": [
    "Do not split cells from one donor across train/validation.",
    "Do not call LD-correlated allele combinations observed haplotypes without phased genotypes.",
    "Do not pool raw junction counts across cohorts/library preparations.",
    "Do not use AlphaGenome outputs to train, calibrate, or define functional priors."
  ]
}
with open(os.path.join(OUT, "preflight_contract.json"), "w") as f: json.dump(contract, f, indent=2)
print(json.dumps({"out": OUT, "resources": contract["resources"]}, indent=2))
