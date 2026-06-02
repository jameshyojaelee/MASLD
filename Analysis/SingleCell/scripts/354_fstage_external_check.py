#!/usr/bin/env python
"""
354_fstage_external_check.py  (S5)

External-cohort F-stage validation for Script 343m augmented inference.

Goal: identify any human scRNA cohort with documented per-donor
F-stage that is NOT in the scVI training set used by 343m. If such a
cohort exists, project its donors through the existing scVI latent,
predict F-stage via kNN-5, and compute quadratic-weighted-kappa (QWK)
against documented stage. Compare to 343m claimed external QWK 0.862.

Strategy:
    1. Re-scan GEO `Sample_characteristics_ch1` for the 7 scRNA datasets
       (GSE136103, GSE174748, GSE185477, GSE189600, GSE244832,
        GSE192740, GSE202379), looking for fibrosis / SAF / METAVIR
       patterns beyond the regex used by 343.
    2. Parse the Wang GSE244832 paper supplementary tables (PDF/XLSX)
       at data/external/wang_gse244832_supp/ if present.
    3. Parse Liver_Atlas (Guilliams GSE192740) supplementary tables at
       data/external/guilliams_liver_atlas/ if present.
    4. Whatever donors come out with documented F-stage AND are not in
       the 343m training (`F_stage_documented` non-NA in
       donor_metadata_extended.tsv -> Andrews-only currently),
       project via scVI latent and predict.
    5. Compute QWK and write fstage_external_report.tsv.

Output: Analysis/SingleCell/results_gpu_v2/fstage_external/
"""

from __future__ import annotations
import os
import re
import sys
import json
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))

EXT_META = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/mcp/inputs/donor_metadata.tsv"
EXTENDED = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/mcp/inputs/donor_metadata_extended.tsv"
OUT_ROOT = Path(os.environ.get("S5_OUT_ROOT", str(PROJECT_ROOT)))
OUT_DIR  = OUT_ROOT / "Analysis/SingleCell/results_gpu_v2/fstage_external"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# GEO matrices scanned by 343; we re-scan more broadly here
GEO_DIRS = {
    "GSE136103": PROJECT_ROOT / "data/raw/GSE136103/geo_metadata",
    "GSE174748": PROJECT_ROOT / "data/raw/GSE174748/geo_metadata",
    "GSE185477": PROJECT_ROOT / "data/raw/GSE185477/geo_metadata",
    "GSE189600": PROJECT_ROOT / "data/raw/GSE189600/geo_metadata",
    "GSE244832": PROJECT_ROOT / "data/raw/GSE244832/geo_metadata",
    "GSE192740": PROJECT_ROOT / "data/raw/GSE192740/geo_metadata",
    "GSE202379": PROJECT_ROOT / "data/raw/GSE202379/geo_metadata",
}

# Broader regex set
FSTAGE_RES = [
    re.compile(r"F\s*([0-4])", re.I),
    re.compile(r"fibrosis\s*stage[:\s]+([0-4])", re.I),
    re.compile(r"METAVIR\s*F\s*([0-4])", re.I),
    re.compile(r"NAS\s*F\s*([0-4])", re.I),
    re.compile(r"Brunt\s*F\s*([0-4])", re.I),
]


def scan_geo_text() -> pd.DataFrame:
    rows = []
    for dataset, gdir in GEO_DIRS.items():
        if not gdir.exists():
            continue
        for txt in gdir.rglob("*"):
            if not txt.is_file():
                continue
            try:
                content = txt.read_text(errors="ignore")[:5_000_000]
            except Exception:
                continue
            for line in content.splitlines():
                if not re.search(r"GSM|donor|sample", line, re.I):
                    continue
                for rgx in FSTAGE_RES:
                    m = rgx.search(line)
                    if m:
                        rows.append({"dataset": dataset,
                                     "line": line[:200],
                                     "F_stage": int(m.group(1)),
                                     "regex": rgx.pattern})
                        break
    return pd.DataFrame(rows)


def main():
    # Load existing donor metadata
    if EXTENDED.exists():
        meta = pd.read_csv(EXTENDED, sep="\t")
    else:
        meta = pd.read_csv(EXT_META, sep="\t")
    print(f"[load] donor metadata: {meta.shape}", flush=True)

    # Document which donors currently have F_stage_documented (Andrews-only
    # expected per CLAUDE.md)
    docs = meta[meta.get("F_stage_documented").notna()] if "F_stage_documented" in meta.columns else pd.DataFrame()
    print(f"[meta] F_stage_documented = {len(docs)} donors", flush=True)

    # Scan GEO text for new F-stage hits
    geo_hits = scan_geo_text()
    geo_hits.to_csv(OUT_DIR / "geo_fstage_scan.tsv", sep="\t", index=False)
    print(f"[scan] geo hits = {len(geo_hits)}", flush=True)

    # Look for supplementary-table dumps
    supp_paths = [
        PROJECT_ROOT / "data/external/wang_gse244832_supp",
        PROJECT_ROOT / "data/external/guilliams_liver_atlas",
        PROJECT_ROOT / "data/external/andrews_2022_supp",
    ]
    supp_present = {p.name: p.exists() for p in supp_paths}
    (OUT_DIR / "supp_availability.json").write_text(
        json.dumps(supp_present, indent=2))
    print(f"[supp] {supp_present}", flush=True)

    # If no NEW documented F-stage donors -> external check is blocked
    # Otherwise, we'd project here. We document the blocker explicitly.
    held_out = pd.DataFrame()
    if "F_stage_documented" in meta.columns and "dataset" in meta.columns:
        # Per CLAUDE.md, currently only Andrews (GSE202379) has documented F.
        non_andrews_docs = meta[
            meta["F_stage_documented"].notna() &
            (meta["dataset"] != "GSE202379")
        ]
        held_out = non_andrews_docs.copy()

    held_out.to_csv(OUT_DIR / "candidate_external_donors.tsv",
                    sep="\t", index=False)
    print(f"[held-out] candidate external donors = {len(held_out)}",
          flush=True)

    # Project + kNN-5 if scvi latent is on disk and we have held-out donors
    scvi_latent = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/scvi/latent.csv.gz"
    train_pred  = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_fstage_scanvi_augmented_predicted.tsv"

    report = {
        "n_donors_F_stage_documented": int(len(docs)),
        "n_external_candidates": int(len(held_out)),
        "geo_scan_hits": int(len(geo_hits)),
        "supp_availability": supp_present,
        "scvi_latent_present": scvi_latent.exists(),
        "train_pred_present": train_pred.exists(),
        "note": ("If n_external_candidates == 0, no documented F-stage cohort "
                 "exists outside the 343m training set. The claimed external "
                 "QWK 0.862 in 343m is therefore an internal LOOCV metric, "
                 "not an external-cohort metric. Genuine external validation "
                 "requires harvesting Wang GSE244832 or Liver_Atlas "
                 "supplementary tables (currently missing on disk)."),
    }
    (OUT_DIR / "external_check_report.json").write_text(
        json.dumps(report, indent=2))
    print(json.dumps(report, indent=2), flush=True)

    # If we have both latent + held-out donors, run kNN-5 QWK.
    if scvi_latent.exists() and len(held_out) > 0 and "F_stage_augmented" in meta.columns:
        try:
            from sklearn.neighbors import KNeighborsClassifier
            from sklearn.metrics import cohen_kappa_score
            L = pd.read_csv(scvi_latent, index_col=0)
            train = meta[meta["F_stage_documented"].notna() &
                         (meta["dataset"] == "GSE202379")]
            X_tr = L.loc[L.index.intersection(train["sample"])].values
            y_tr = train.set_index("sample").loc[
                L.index.intersection(train["sample"]), "F_stage_documented"].astype(int).values
            X_te = L.loc[L.index.intersection(held_out["sample"])].values
            y_te = held_out.set_index("sample").loc[
                L.index.intersection(held_out["sample"]), "F_stage_documented"].astype(int).values
            if len(X_tr) > 5 and len(X_te) > 1:
                knn = KNeighborsClassifier(n_neighbors=5)
                knn.fit(X_tr, y_tr)
                pred = knn.predict(X_te)
                qwk = cohen_kappa_score(y_te, pred, weights="quadratic")
                report["external_qwk_knn5"] = float(qwk)
                report["external_n"] = int(len(y_te))
                (OUT_DIR / "external_check_report.json").write_text(
                    json.dumps(report, indent=2))
                print(f"[external QWK kNN-5] = {qwk:.3f} on n={len(y_te)}",
                      flush=True)
        except Exception as e:
            print(f"[external projection failed] {e}", flush=True)


if __name__ == "__main__":
    main()
