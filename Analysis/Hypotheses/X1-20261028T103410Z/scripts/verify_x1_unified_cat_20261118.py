#!/usr/bin/env python3
"""Verify X1 rerun with Unified Gene Catalog 2026-11-18 and attach refreshed Catalog annotations.

Checks that running x1_xenium_v2_20261118.py with X1_CATALOG_PATH pointed to
gene-catalog-unified-20261118 produces ZERO carrier, effect-size, p-value, or CI changes
relative to out_20261118, and writes Catalog-annotated companion tables for results_5d,
results_5e, and results_calls_supplementary.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OLD_DIR = ROOT / "Analysis/Hypotheses/X1-20261028T103410Z/out_20261118"
NEW_DIR = Path(sys.argv[1])
CAT_PATH = Path(sys.argv[2])

cat = pd.read_csv(
    CAT_PATH, sep="\t",
    usecols=[
        "symbol", "ensembl_id",
        "baseline_primary_evidence_class_20261106",
        "primary_evidence_class",
        "candidate_class",
        "genetic_mechanism_candidate",
        "disease_liability_candidate",
        "coding_fine_mapped",
        "coding_tier2_enzyme_only",
        "state_partition",
        "state_partition_c3",
        "state_hep_ambient_corrected",
    ],
).drop_duplicates("symbol").set_index("symbol")

audit_rows = []
for fname in ("results_5d.tsv", "results_5e.tsv", "results_calls_supplementary.tsv"):
    old_df = pd.read_csv(OLD_DIR / fname, sep="\t")
    new_df = pd.read_csv(NEW_DIR / fname, sep="\t")
    assert list(old_df.columns) == list(new_df.columns), f"Columns differ in {fname}"
    assert len(old_df) == len(new_df), f"Row count differs in {fname}: {len(old_df)} vs {len(new_df)}"
    assert (old_df["gene"].values == new_df["gene"].values).all(), f"Gene order differs in {fname}"
    assert (old_df["carrier"].values == new_df["carrier"].values).all(), f"Carrier differs in {fname}"

    max_num_diff = 0.0
    str_diffs = 0
    for col in old_df.columns:
        if pd.api.types.is_numeric_dtype(old_df[col]) and pd.api.types.is_numeric_dtype(new_df[col]):
            a = old_df[col].to_numpy(dtype=float)
            b = new_df[col].to_numpy(dtype=float)
            na_same = np.isnan(a) == np.isnan(b)
            assert na_same.all(), f"NaN mismatch in {fname}:{col}"
            inf_same = (np.isposinf(a) == np.isposinf(b)) & (np.isneginf(a) == np.isneginf(b))
            assert inf_same.all(), f"Inf mismatch in {fname}:{col}"
            mask = np.isfinite(a) & np.isfinite(b)
            if mask.any():
                d = float(np.max(np.abs(a[mask] - b[mask])))
                if d > max_num_diff:
                    max_num_diff = d
                assert d < 1e-12, f"Numeric diff {d} in {fname}:{col}"
        else:
            sa = old_df[col].fillna("__NA__").astype(str)
            sb = new_df[col].fillna("__NA__").astype(str)
            if not (sa == sb).all():
                str_diffs += int((sa != sb).sum())
                raise AssertionError(f"String mismatch in {fname}:{col}")

    ann = new_df.merge(cat.reset_index(), left_on="gene", right_on="symbol", how="left").drop(columns=["symbol"])
    out_ann = NEW_DIR / fname.replace(".tsv", "_catalog_annotated.tsv")
    ann.to_csv(out_ann, sep="\t", index=False)

    n_reclass = int((ann["baseline_primary_evidence_class_20261106"] != ann["candidate_class"]).sum())
    audit_rows.append(dict(
        file=fname, n_genes=len(new_df),
        max_numeric_diff_vs_out_20261118=max_num_diff,
        string_diffs_vs_out_20261118=str_diffs,
        n_matched_in_unified_catalog=int(ann["candidate_class"].notna().sum()),
        n_diff_from_baseline_20261106_incl_rename=n_reclass,
        candidate_class_counts=";".join(f"{k}:{v}" for k, v in ann["candidate_class"].value_counts(dropna=False).items()),
    ))

audit = pd.DataFrame(audit_rows)
audit.to_csv(NEW_DIR / "catalog_refresh_audit.tsv", sep="\t", index=False)
print(audit.to_string(index=False))
print("X1_UNIFIED_CATALOG_VERIFY_PASS")
