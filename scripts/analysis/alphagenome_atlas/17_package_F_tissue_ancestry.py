#!/usr/bin/env python3
"""Step 17 (Package F): tissue-panel sensitivity and ancestry robustness.

Tissue: the fixed panel (tissue_panel.json) on identical scorer definitions; tracks summarised within
group by median before groups are compared (no max-across-tissues advantage). HepG2 kept separate.
Ancestry: the same static annotation applied to each ancestry-specific direct-trait credible set
(uniform35_v4: EUR, AMR, AFR studies of the same trait); differences decomposed into posterior
(which variants), coverage, and contributing-study differences. The same variant receives the same
static prediction whatever the ancestry.

Outputs (tables/): tissue_panel_signal_summary.tsv, ancestry_comparison.tsv
"""

from __future__ import annotations

import csv
import json
import math
from collections import defaultdict

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

import lib_atlas as la

TABLES = la.out_root() / "tables"
PANEL = json.load((la.SCRIPT_DIR / "tissue_panel.json").open())
SCORERS = ["RNA_SEQ", "ATAC", "DNASE", "CHIP_HISTONE"]


def main() -> None:
    signals = {s["signal_uid"]: s for s in la.read_tsv(TABLES / "eligible_signals.tsv")}
    weights = defaultdict(dict)
    with la.open_text(TABLES / "signal_variant_weights.tsv.gz") as h:
        for r in csv.DictReader(h, delimiter="\t"):
            if r["mapping_status"] == "mapped" and r["in_query_set"] == "True":
                weights[r["signal_uid"]][r["variant_uid"]] = float(r["weight"])
    rows = []
    for scorer in SCORERS:
        path = TABLES / "prediction_records" / f"{scorer}.parquet"
        if not path.exists():
            continue
        t = pq.read_table(path).to_pandas()
        groups = sorted({c.split("|")[1] for c in t.columns if c.startswith("q|")})
        # per row: median signed quantile per group
        med = {g: t[[c for c in t.columns if c.startswith(f"q|{g}|")]].median(axis=1).values for g in groups}
        t_idx = t.set_index(["variant_uid", "gene_id"]).index if "gene_id" in t else t.set_index("variant_uid").index
        for sig, w in weights.items():
            s = signals[sig]
            target = s["ensembl"]
            if "gene_id" in t and target:
                mask = (t["variant_uid"].isin(w)) & (t["gene_id"] == target)
            else:
                mask = t["variant_uid"].isin(w) & (t["gene_id"] == "" if "gene_id" in t else True)
            sub = t[mask]
            if len(sub) == 0:
                continue
            wv = np.array([w[v] for v in sub["variant_uid"]])
            cov = wv.sum()
            row = {"signal_uid": sig, "universe": s["universe"], "scorer": scorer, "gene_scope": "target_gene" if ("gene_id" in t and target) else "region",
                   "coverage": cov, "n_variants": len(sub)}
            for g in groups:
                vals = med[g][mask.values]
                ok = ~np.isnan(vals)
                row[f"{g}_A_quantile_over_C"] = float((wv[ok] * np.abs(vals[ok])).sum() / wv[ok].sum()) if ok.any() else math.nan
                row[f"{g}_S_quantile_over_C"] = float((wv[ok] * vals[ok]).sum() / wv[ok].sum()) if ok.any() else math.nan
            liver = row.get("primary_liver_A_quantile_over_C", math.nan)
            others = {g: row.get(f"{g}_A_quantile_over_C", math.nan) for g in groups if g not in ("primary_liver", "hepatocyte", "HepG2")}
            row["liver_specific"] = (not math.isnan(liver)) and liver >= 0.9 and all((math.isnan(v) or v < 0.5) for v in others.values())
            row["shared_across_panel"] = (not math.isnan(liver)) and liver >= 0.9 and any((not math.isnan(v)) and v >= 0.9 for v in others.values())
            row["note"] = "predicted activity outside liver is not evidence of extrahepatic mediation"
            rows.append(row)
    la.write_tsv_once(TABLES / "tissue_panel_signal_summary.tsv", rows, sorted({k for r in rows for k in r}))

    # ancestry: direct-trait credible sets per trait across ancestries (Universe B) — same block, different ancestry
    b = [s for s in signals.values() if s["universe"] == "B_direct"]
    ps = pd.read_csv(TABLES / "posterior_summaries.tsv.gz", sep="\t")
    rna = ps[(ps["scorer"] == "RNA_SEQ") & (ps["gene_scope"].isin(["target_gene", "region"]))].drop_duplicates("signal_uid").set_index("signal_uid")
    atac = ps[(ps["scorer"] == "ATAC")].drop_duplicates("signal_uid").set_index("signal_uid")
    by_block = defaultdict(list)
    for s in b:
        by_block[(s["analysis_block"], s["trait"])].append(s)
    anc_rows = []
    for (block, trait), group in by_block.items():
        ancs = {s["ancestry"] for s in group}
        if len(ancs) < 2:
            continue
        for s in group:
            sig = s["signal_uid"]
            w = weights[sig]
            top = max(w, key=w.get) if w else ""
            anc_rows.append({"analysis_block": block, "trait": trait, "signal_uid": sig, "study": s["gwas_name"], "ancestry": s["ancestry"],
                             "n_cs_variants": s["n_cs_variants"], "cs_coverage": s["cs_coverage"], "top_variant": top, "top_pip": w.get(top, math.nan),
                             "queried_mass": sum(w.values()),
                             "atac_C": atac.at[sig, "C_liver"] if sig in atac.index else math.nan, "atac_A_quantile": atac.at[sig, "A_quantile_weighted"] if sig in atac.index else math.nan,
                             "shared_top_variant_with_other_ancestries": any(top in weights[o["signal_uid"]] for o in group if o["signal_uid"] != sig),
                             "decomposition": "same variant, same static prediction; differences arise from posterior weights, coverage and contributing studies"})
    la.write_tsv_once(TABLES / "ancestry_comparison.tsv", anc_rows, list(anc_rows[0].keys()) if anc_rows else ["analysis_block"])
    la.log(f"Package F: {len(rows)} tissue rows, {len(anc_rows)} ancestry rows")


if __name__ == "__main__":
    main()
