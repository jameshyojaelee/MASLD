#!/usr/bin/env python
"""Assemble every number the IGFBP7 ambient SENSITIVITY ANALYSIS verdict rests on.

SCOPE (binding): reads this candidate's own outputs plus the read-only frozen
registry. Writes only into this candidate. Nothing frozen is modified.
"""
from __future__ import annotations

import hashlib
import os
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(os.environ["MASLD_PROJECT_ROOT"])
CAND = Path(os.environ["CAND_ROOT"])
RES = CAND / "results"
PROGRAM_ROOT = (
    ROOT
    / "Analysis/Multimodal_Program_Projection/candidates"
    / "program-context-v2-candidate-2026-08-07/hotspot"
)

FOCUS = {
    8: "Stromal ECM (IGFBP7) - UNDER TEST",
    20: "Ductular injury (BICC1) - paired comparator",
    2: "Sinusoidal endothelial (STAB2) - known contamination positive control",
    1: "Leukocyte / immune (DOCK2) - known contamination positive control",
    12: "Secretory plasma protein (ALB) - hepatocyte-intrinsic negative control",
    21: "Xenobiotic / drug metab (CYP) - hepatocyte-intrinsic negative control",
    16: "PPARa lipid metabolism - hepatocyte-intrinsic negative control",
}


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main() -> None:
    eff = pd.read_csv(RES / "03_program_stage_effects.tsv", sep="\t")
    burden = pd.read_csv(RES / "04_program_ambient_burden.tsv", sep="\t")
    cal = pd.read_csv(RES / "05_own_source_calibration.tsv", sep="\t")

    m = eff.merge(
        burden[
            ["module", "l1_weighted_ambient_fraction_hep",
             "program_count_ambient_fraction_hep", "frac_members_ambient_gt_0.8"]
        ],
        on="module", how="left",
    ).merge(
        cal[["module", "source_lineage", "spearman_vs_own_source",
             "spearman_vs_fibroblast", "spearman_vs_non_hepatocyte"]],
        on="module", how="left",
    )
    m["stage_beta_pct_retained"] = 100 * m["corrected_beta"] / m["refit_stored_beta"]
    m = m.sort_values("l1_weighted_ambient_fraction_hep", ascending=False)
    m.to_csv(RES / "06_master_summary.tsv", sep="\t", index=False)

    print("=" * 118)
    print("ALL 30 FROZEN HEPATOCYTE PROGRAMS: ambient burden vs stage-effect survival")
    print("=" * 118)
    cols = [
        "module", "module_name", "robust_display",
        "l1_weighted_ambient_fraction_hep", "donor_pearson_stored_vs_corrected",
        "refit_stored_beta", "corrected_beta", "stage_beta_pct_retained",
        "corrected_hc3_pvalue", "corrected_qvalue", "direction_preserved",
    ]
    with pd.option_context("display.width", 250, "display.max_columns", 40):
        print(m[cols].to_string(index=False, float_format=lambda v: f"{v:.4g}"))

    print()
    print("=" * 118)
    print("FOCUS PROGRAMS")
    print("=" * 118)
    for mod, label in FOCUS.items():
        r = m[m["module"] == mod]
        if r.empty:
            continue
        r = r.iloc[0]
        print(f"\n--- hep module {mod}: {label}")
        print(f"    frozen registry q                 {r['stored_registry_qvalue']:.4g}"
              f"   robust_display={r['robust_display']}")
        print(f"    L1-weighted ambient fraction      {r['l1_weighted_ambient_fraction_hep']:.4f}"
              f"   ({100*r['frac_members_ambient_gt_0.8']:.0f}% of members >0.8)")
        print(f"    own-source lineage rho            {r['spearman_vs_own_source']}"
              f"   (source={r['source_lineage']})")
        print(f"    donor score r / rho stored vs corrected   "
              f"{r['donor_pearson_stored_vs_corrected']:.4f} / "
              f"{r['donor_spearman_stored_vs_corrected']:.4f}")
        print(f"    stage beta   {r['refit_stored_beta']:+.4f} -> {r['corrected_beta']:+.4f}"
              f"   ({r['stage_beta_pct_retained']:.1f}% retained, direction preserved="
              f"{r['direction_preserved']})")
        print(f"    stage SE     {r['refit_stored_se']:.4f} -> {r['corrected_se']:.4f}"
              f"   HC3 SE {r['refit_stored_hc3_se']:.4f} -> {r['corrected_hc3_se']:.4f}")
        print(f"    p            {r['refit_stored_pvalue']:.4g} -> {r['corrected_pvalue']:.4g}"
              f"   HC3 p {r['refit_stored_hc3_pvalue']:.4g} -> {r['corrected_hc3_pvalue']:.4g}")
        print(f"    q (BH over 117)   {r['stored_registry_qvalue']:.4g} -> "
              f"{r['corrected_qvalue']:.4g}   HC3 q -> {r['corrected_hc3_qvalue']:.4g}")

    print()
    print("=" * 118)
    print("CONTROL SEPARATION CHECK (does correction move everything equally?)")
    print("=" * 118)
    known = m[m["module"].isin([1, 2])]
    intrinsic = m[m["module"].isin([12, 21, 16])]
    test = m[m["module"] == 8]
    comp = m[m["module"] == 20]
    for name, grp in [("known contamination (1,2)", known),
                      ("IGFBP7 (8)", test),
                      ("BICC1 (20)", comp),
                      ("hepatocyte-intrinsic (12,21,16)", intrinsic),
                      ("all 30", m)]:
        print(f"  {name:<34} mean ambient={grp['l1_weighted_ambient_fraction_hep'].mean():.4f}"
              f"  mean %beta retained={grp['stage_beta_pct_retained'].mean():8.2f}"
              f"  mean donor r={grp['donor_pearson_stored_vs_corrected'].mean():.4f}")

    # output checksum manifest
    rows = []
    for p in sorted(RES.glob("*")):
        if p.is_file():
            rows.append({"role": "output", "path": str(p), "bytes": p.stat().st_size,
                         "sha256": sha256(p)})
    for p in [PROGRAM_ROOT / "program_registry_v2.tsv",
              PROGRAM_ROOT / "program_membership_v2.tsv"]:
        rows.append({"role": "frozen_input_readonly", "path": str(p),
                     "bytes": p.stat().st_size, "sha256": sha256(p)})
    pd.DataFrame(rows).to_csv(RES / "06_output_checksums.tsv", sep="\t", index=False)
    print("\n[assemble] wrote 06_master_summary.tsv and 06_output_checksums.tsv", flush=True)


if __name__ == "__main__":
    main()
