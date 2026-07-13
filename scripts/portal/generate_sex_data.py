#!/usr/bin/env python3
"""Sex layer web data (data contract §5.5).

Source (PINNED): sex_v3/lvqw_staging/sex_deg_classification_v3.csv, `sex_class`
col = exactly 7 Male_biased + 1 Female_biased (the LVQW-collapse headline 8).
Do NOT use the top-level v3 (11F/4M) or the atlas sex_class (6).

Framing: LVQW-fixed (C2) sex layer; the female-progressor narrative was retired
2026-06-08. Sex-biased framing only.

Outputs (-> --output-dir):
  sex_summary.json
  sex_gene_classification.parquet   (8 rows, sorted by symbol)

Run:
  micromamba run -n spatial python scripts/portal/generate_sex_data.py
"""

from __future__ import annotations

import argparse

import pandas as pd

from _portal_io import dump_json, project_root, r, write_parquet

NOTE = ("LVQW-fixed (C2) sex layer; female-progressor narrative retired "
        "2026-06-08. Sex-biased framing only.")


def main():
    root = project_root()
    ap = argparse.ArgumentParser()
    ap.add_argument("--output-dir", default=str(root / "masld-atlas-v2/public/data"))
    args = ap.parse_args()
    out = __import__("pathlib").Path(args.output_dir)

    src = (root / "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/"
           "integration/sex_v3/lvqw_staging/sex_deg_classification_v3.csv")
    print(f"Loading {src}")
    df = pd.read_csv(src)
    df = df[df["sex_class"].isin(["Male_biased", "Female_biased"])].copy()

    n_male = int((df["sex_class"] == "Male_biased").sum())
    n_female = int((df["sex_class"] == "Female_biased").sum())
    assert len(df) == 8 and n_male == 7 and n_female == 1, (
        f"expected 8 dimorphic genes (7M/1F), got {len(df)} ({n_male}M/{n_female}F)")

    rows = []
    for _, x in df.iterrows():
        rows.append({
            "symbol": str(x["gene_symbol"]),
            "ensembl": str(x["ensembl_base"]),
            "chr": str(x["chr"]),
            "biotype": str(x["gene_biotype"]),
            "sex_class": str(x["sex_class"]),
            "logfc_female": r(x["logFC_F"], 4),
            "logfc_male": r(x["logFC_M"], 4),
            "interaction_padj": r(x["interaction_padj"], 6),
            "beta_female": r(x["beta_F"], 4),
            "beta_male": r(x["beta_M"], 4),
        })
    rows.sort(key=lambda d: d["symbol"])

    summary = {
        "n_dimorphic": len(rows),
        "n_female_biased": n_female,
        "n_male_biased": n_male,
        "genes": rows,
        "note": NOTE,
    }
    kb = dump_json(summary, out / "sex_summary.json")
    print(f"  -> sex_summary.json: {kb:.1f} KB ({n_male}M/{n_female}F)")

    pq = pd.DataFrame(rows)[[
        "symbol", "ensembl", "chr", "biotype", "sex_class",
        "logfc_female", "logfc_male", "interaction_padj",
        "beta_female", "beta_male"]]
    write_parquet(pq, out / "sex_gene_classification.parquet")
    print(f"  -> sex_gene_classification.parquet: {len(pq)} rows")


if __name__ == "__main__":
    main()
