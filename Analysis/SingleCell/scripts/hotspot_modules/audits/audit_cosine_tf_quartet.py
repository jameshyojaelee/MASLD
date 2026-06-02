"""Audit the cosine=1.0 artifact in novelty_matches.tsv and check the
defensibility of the 4 TF labels (hep__10/RXRA, fib__13/RORA, endo__22/THRB,
mac__30/Hedgehog).

Findings to produce:
  - Distribution of best_match_cosine across 202 modules
  - Root cause of cosine=1.0 (intersection-based cosine with constant-1 ref)
  - Per-quartet: TF-in-top-50, clean Jaccard@50, hypergeometric p, verdict
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import hypergeom

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
RESULTS = ROOT / "Analysis/SingleCell/results_gpu_v2/hotspot_modules"
REF = ROOT / "data/reference_signatures"
OUT = ROOT / "Analysis/SingleCell/scripts/hotspot_modules/audits"

# ---- Q1: cosine distribution ----
nm = pd.read_csv(RESULTS / "novelty_matches.tsv", sep="\t")
print(f"Total modules: {len(nm)}")
c = nm["best_match_cosine"]
print(f"cosine == 1.0 exact: {(c == 1.0).sum()}")
print(f"cosine >= 0.999     : {(c >= 0.999).sum()}")
print(f"cosine >= 0.95      : {(c >= 0.95).sum()}")
print(f"cosine >= 0.50      : {(c >= 0.50).sum()}")
print(f"cosine <  0.35      : {(c < 0.35).sum()}")

# ---- Q1: how often is best-match panel a membership panel (constant 1.0)? ----
print("\nbest_match_panel distribution among cosine>=0.99 rows:")
print(nm.loc[c >= 0.99, "best_match_panel"].value_counts())

# ---- Q2: load quartet modules and SCENIC/Hallmark refs ----
quartet = [
    ("hepatocytes", 10, "scenic_hep", "RXRA"),
    ("fibroblasts", 13, "scenic_hep", "RORA"),
    ("endothelial_cells", 22, "scenic_hep", "THRB"),
    ("macrophages", 30, "hallmark", "HALLMARK_HEDGEHOG_SIGNALING"),
]

scenic = pd.read_csv(REF / "scenic_hep_regulons.tsv", sep="\t")
hallmark = pd.read_csv(REF / "hallmark_v2025_1.tsv", sep="\t")
ref_all = pd.concat([scenic, hallmark], ignore_index=True)

# Universe = union of all module genes + all reference genes (proxy for tested universe)
module_universe = set()
quartet_rows = []
for ct, mod, panel, prog in quartet:
    mg = pd.read_csv(RESULTS / ct / "module_genes.tsv", sep="\t")
    module_universe |= set(mg["gene"].unique())
    sub = mg.loc[mg["module"] == mod].sort_values("weight", ascending=False)
    top50 = sub.head(50)["gene"].tolist()
    all_module_genes = set(sub["gene"])

    if panel == "scenic_hep":
        ref_genes = set(scenic.loc[scenic["program"] == prog, "gene"])
    else:
        ref_genes = set(hallmark.loc[hallmark["program"] == prog, "gene"])

    tf_name = prog if panel == "scenic_hep" else None
    tf_in_top50 = (tf_name in top50) if tf_name else None
    tf_in_module = (tf_name in all_module_genes) if tf_name else None

    inter_top50 = ref_genes & set(top50)
    union_top50 = ref_genes | set(top50)
    jac50 = len(inter_top50) / max(len(union_top50), 1)

    inter_all = ref_genes & all_module_genes
    jac_all = len(inter_all) / max(len(ref_genes | all_module_genes), 1)

    # Clean cosine with union zero-padding (binary ref, weighted module)
    union_idx = sorted(all_module_genes | ref_genes)
    w_module = mg.loc[mg["module"] == mod].set_index("gene")["weight"]
    u = np.array([w_module.get(g, 0.0) for g in union_idx])
    v = np.array([1.0 if g in ref_genes else 0.0 for g in union_idx])
    nu, nv = np.linalg.norm(u), np.linalg.norm(v)
    clean_cos = float(np.dot(u, v) / (nu * nv)) if nu * nv > 0 else 0.0

    quartet_rows.append({
        "module": f"{ct}__{mod}",
        "label": prog,
        "n_module_genes": len(all_module_genes),
        "n_ref_genes": len(ref_genes),
        "tf_in_top50": tf_in_top50,
        "tf_in_module": tf_in_module,
        "intersect_top50": ",".join(sorted(inter_top50)) or "(none)",
        "intersect_all": ",".join(sorted(inter_all)) or "(none)",
        "jaccard_top50_clean": round(jac50, 4),
        "jaccard_all_clean": round(jac_all, 4),
        "clean_cosine_zeropad": round(clean_cos, 4),
    })

# Hypergeometric p-value: universe = module_universe ∪ ref_genes_all
universe = module_universe | set(ref_all["gene"])
N = len(universe)
print(f"\nGene universe size for hypergeometric: {N}")

for row in quartet_rows:
    ct, mod = row["module"].split("__")
    mod = int(mod)
    panel = "scenic_hep" if row["label"].upper() == row["label"] and not row["label"].startswith("HALLMARK") else "hallmark"
    if row["label"].startswith("HALLMARK"):
        ref_genes = set(hallmark.loc[hallmark["program"] == row["label"], "gene"])
    else:
        ref_genes = set(scenic.loc[scenic["program"] == row["label"], "gene"])
    K = len(ref_genes & universe)
    n = 50  # top-50 draw
    mg = pd.read_csv(RESULTS / ct / "module_genes.tsv", sep="\t")
    top50 = set(mg.loc[mg["module"] == mod].sort_values("weight", ascending=False).head(50)["gene"])
    k = len(top50 & ref_genes)
    # P(X >= k) under hypergeometric(N, K, n)
    p = hypergeom.sf(k - 1, N, K, n) if k > 0 else 1.0
    row["hyper_k_intersect"] = k
    row["hyper_K_refsize"] = K
    row["hyper_pval"] = float(f"{p:.3g}")

audit = pd.DataFrame(quartet_rows)
print("\n=== TF quartet audit ===")
print(audit.to_string(index=False))
audit.to_csv(OUT / "tf_quartet_audit.tsv", sep="\t", index=False)

# Save cosine distribution summary
dist = pd.DataFrame({
    "bucket": ["==1.0", ">=0.999", ">=0.95", ">=0.50", "<0.35"],
    "n_modules": [
        int((c == 1.0).sum()),
        int((c >= 0.999).sum()),
        int((c >= 0.95).sum()),
        int((c >= 0.50).sum()),
        int((c < 0.35).sum()),
    ],
})
dist.to_csv(OUT / "cosine_distribution.tsv", sep="\t", index=False)
print("\n=== Cosine distribution ===")
print(dist.to_string(index=False))
print(f"\nWrote {OUT}/tf_quartet_audit.tsv and {OUT}/cosine_distribution.tsv")
