"""Read full-run + LOO module_genes.tsv files, compute top-50 Jaccard stability.

Run this after ALL 49 LOO tasks (503_loo_stability.sbatch) have completed.

Output: results_gpu_v2/hotspot_modules/loo_stability.tsv
  Columns: cell_type, module, stability_score, wang_loo_jaccard (NaN — no Wang in this
           atlas; GSE244832 is the large Wang-equivalent dataset), loo_<dataset>, ...,
           stability_fail (stability_score < 0.5)
"""
from __future__ import annotations
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from hotspot_io import RESULTS, RUN_ORDER

# Actual dataset names from donor_metadata.tsv `dataset` column (enumerated 2026-05-13)
DATASETS = [
    "GSE136103",
    "GSE174748",
    "GSE185477",
    "GSE189600",
    "GSE202379",
    "GSE244832",
    "Liver_Atlas",
]

# Post-2026-05-22 protocol-contamination remediation: the re-discovered non-hep
# cell types (mac/fib/chol/tcells) were fit with GSE136103 + Liver_Atlas excluded,
# and their clean LOO re-runs (503_loo_clean.sbatch) hold out only these 5. Reading
# the stale/contaminated holdouts for those CTs would mis-pair pre- vs post-
# remediation modules, so restrict them to the clean set. Hepatocytes were NOT
# re-discovered and legitimately use the Liver_Atlas hepatocyte fraction
# (GSE136103 contributes ~0 hep cells), so they keep the full dataset list.
CLEAN_DATASETS = [
    "GSE174748",
    "GSE185477",
    "GSE189600",
    "GSE202379",
    "GSE244832",
]


def datasets_for(ct: str) -> list[str]:
    return DATASETS if ct == "hepatocytes" else CLEAN_DATASETS


def top50(df: pd.DataFrame, module: int) -> set[str]:
    """Return the top-50 genes (by weight) for a given module."""
    return set(
        df[df["module"] == module]
        .sort_values("weight", ascending=False)
        .head(50)["gene"]
    )


def main() -> None:
    rows = []
    for ct in RUN_ORDER:
        full_path = RESULTS / ct / "module_genes.tsv"
        if not full_path.exists():
            print(f"[SKIP] {ct}: full-run module_genes.tsv not found — skipping.")
            continue
        full = pd.read_csv(full_path, sep="\t")

        for mod in sorted(full["module"].unique()):
            top_full = top50(full, mod)
            if not top_full:
                continue

            per_loo: dict[str, float] = {}
            for ds in datasets_for(ct):
                loo_path = RESULTS / ct / "loo" / ds / "module_genes.tsv"
                if not loo_path.exists():
                    per_loo[ds] = float("nan")
                    continue
                loo = pd.read_csv(loo_path, sep="\t")
                # Best-match Jaccard across all LOO modules (modules may be renumbered)
                best = 0.0
                for loo_mod in loo["module"].unique():
                    t = top50(loo, loo_mod)
                    if not t:
                        continue
                    j = len(top_full & t) / max(len(top_full | t), 1)
                    if j > best:
                        best = j
                per_loo[ds] = best

            stability = pd.Series(per_loo).mean(skipna=True)

            rows.append(
                {
                    "cell_type": ct,
                    "module": mod,
                    "stability_score": round(stability, 4),
                    # GSE244832 is the largest single dataset (Wang-equivalent);
                    # report separately to flag sensitivity to this dominant cohort
                    "gse244832_loo_jaccard": per_loo.get("GSE244832", float("nan")),
                    **{f"loo_{ds}": v for ds, v in per_loo.items()},
                    "stability_fail": stability < 0.5,
                }
            )

    if not rows:
        print("No results found — have all LOO runs completed?")
        sys.exit(1)

    out = pd.DataFrame(rows)
    out_path = RESULTS / "loo_stability.tsv"
    out.to_csv(out_path, sep="\t", index=False)
    n_fail = int(out["stability_fail"].sum())
    print(
        f"Wrote {len(out)} module-stability rows to {out_path}\n"
        f"  stability_fail (<0.50): {n_fail}/{len(out)}"
    )


if __name__ == "__main__":
    main()
