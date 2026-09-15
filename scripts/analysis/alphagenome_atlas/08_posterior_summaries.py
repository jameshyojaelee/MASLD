#!/usr/bin/env python3
"""Step 08: §5 posterior summaries per signal × scorer × (gene) over liver tracks.

Coverage is reported in separate layers:
  C_map      mass of variants with a valid hg38 SNV mapping
  C_query    mass of mapped variants inside the prespecified query set
  C_atlas    mass of queried variants returned by the Atlas for this scorer
  C_gene     (gene scorers) mass of returned variants carrying a row for the signal's colocalized gene
  C_liver    mass of variants with ≥1 primary-liver/hepatocyte track for this scorer (= the summary's coverage)
Signed/magnitude summaries use the liver median across tracks; unsigned scorers get magnitude only.
A variant enters a signal once; the join is asserted 1:1 on (signal, variant, scorer, gene).
"""

from __future__ import annotations

import csv
import math
from collections import defaultdict

import pandas as pd

import lib_atlas as la

P = la.prespec()
ROOT = la.out_root()
TABLES = ROOT / "tables"
CUT = P["coverage"]["completeness_cutoff"]

OUT_COLS = ["signal_uid", "universe", "gwas_name", "trait_class", "gene_target", "scorer", "is_signed", "gene_id", "gene_scope",
            "n_variants", "n_mapped", "n_queried", "n_returned", "n_liver_scored", "C_map", "C_query", "C_atlas", "C_gene", "C_liver",
            "S_liver_median", "A_liver_median", "S_over_C", "A_over_C", "S_quantile_weighted", "A_quantile_weighted",
            "missing_mass", "conditional_flag", "complete_0.95", "top_variant", "top_weight", "top_variant_liver_median",
            "posterior_definition"]


def main() -> None:
    signals = {s["signal_uid"]: s for s in la.read_tsv(TABLES / "eligible_signals.tsv")}
    meta = {r["scorer"]: r for r in la.read_tsv(TABLES / "scorer_metadata.tsv")}
    # signal -> variant weights (all), mapped, queried
    w_all, w_map, w_q = defaultdict(dict), defaultdict(dict), defaultdict(dict)
    with la.open_text(TABLES / "signal_variant_weights.tsv.gz") as handle:
        for r in csv.DictReader(handle, delimiter="\t"):
            sig, w = r["signal_uid"], float(r["weight"])
            key = r["variant_uid"] or r["source_variant_id"]
            if key in w_all[sig]:
                raise la.ContractError(f"duplicate variant within signal {sig}: {key}")
            w_all[sig][key] = w
            if r["mapping_status"] == "mapped":
                w_map[sig][r["variant_uid"]] = w
                if r["in_query_set"] == "True":
                    w_q[sig][r["variant_uid"]] = w
    avail = pd.read_csv(TABLES / "atlas_availability.tsv.gz", sep="\t", index_col="variant_uid")
    liver = pd.read_csv(TABLES / "liver_summaries.tsv.gz", sep="\t", dtype={"gene_id": str, "gene_name": str}, keep_default_na=False)
    liver["gene_id"] = liver["gene_id"].fillna("").astype(str)
    liver_by_scorer = {s: df for s, df in liver.groupby("scorer")}

    rows = []
    for sig, s in signals.items():
        wa, wm, wq = w_all[sig], w_map[sig], w_q[sig]
        total = sum(wa.values())
        c_map, c_query = sum(wm.values()), sum(wq.values())
        target = s["ensembl"]
        for scorer, m in meta.items():
            is_signed = m["is_signed"] == "True"
            returned = {v: w for v, w in wq.items() if v in avail.index and avail.at[v, scorer] > 0}
            c_atlas = sum(returned.values())
            df = liver_by_scorer.get(scorer)
            if df is None:
                continue
            sub = df[df["variant_uid"].isin(returned)]
            gene_scoped = bool(len(sub)) and (sub["gene_id"] != "").any()
            scopes = []
            if gene_scoped:
                if target:
                    scopes.append(("target_gene", sub[sub["gene_id"] == target]))
                # per-gene rows for every returned gene (Package A multiple-candidate-gene flag)
                for gid, g in sub.groupby("gene_id"):
                    if gid and gid != target:
                        scopes.append(("other_gene", g))
            else:
                scopes.append(("region", sub))
            for scope, g in scopes:
                if len(g) == 0:
                    continue
                dup = g.duplicated(subset=["variant_uid", "gene_id", "junction_start", "junction_end"]).any()
                if dup and scorer != "SPLICE_JUNCTIONS":
                    raise la.ContractError(f"duplicate liver summary rows for {sig} {scorer}")
                # junction scorers: several junctions per variant×gene → collapse to the max-|score| junction (magnitude only)
                if scorer == "SPLICE_JUNCTIONS":
                    g = g.reindex(g["liver_median_raw"].abs().groupby(g["variant_uid"]).idxmax().values)
                scores = dict(zip(g["variant_uid"], g["liver_median_raw"]))
                quant = dict(zip(g["variant_uid"], g["liver_median_quantile"]))
                summ = la.posterior_summary(returned, scores, is_signed, CUT)
                qs = la.posterior_summary(returned, quant, is_signed, CUT)
                c_gene = sum(returned[v] for v in scores) if gene_scoped else math.nan
                gid = g["gene_id"].iloc[0] if gene_scoped else ""
                top = max(returned.items(), key=lambda kv: (kv[1], kv[0]))[0] if returned else ""
                rows.append({
                    "signal_uid": sig, "universe": s["universe"], "gwas_name": s["gwas_name"], "trait_class": s["trait_class"],
                    "gene_target": target, "scorer": scorer, "is_signed": is_signed, "gene_id": gid, "gene_scope": scope,
                    "n_variants": len(wa), "n_mapped": len(wm), "n_queried": len(wq), "n_returned": len(returned), "n_liver_scored": len(scores),
                    "C_map": c_map, "C_query": c_query, "C_atlas": c_atlas, "C_gene": c_gene, "C_liver": summ["coverage"],
                    "S_liver_median": summ["signed"], "A_liver_median": summ["magnitude"], "S_over_C": summ["signed_over_coverage"],
                    "A_over_C": summ["magnitude_over_coverage"], "S_quantile_weighted": qs["signed"], "A_quantile_weighted": qs["magnitude"],
                    "missing_mass": total - summ["coverage"], "conditional_flag": summ["coverage"] < total - 1e-9,
                    "complete_0.95": summ["coverage"] >= CUT, "top_variant": top, "top_weight": returned.get(top, math.nan),
                    "top_variant_liver_median": scores.get(top, math.nan), "posterior_definition": s["posterior_definition"],
                })
    n = la.write_tsv_once(TABLES / "posterior_summaries.tsv.gz", rows, OUT_COLS)
    la.log(f"wrote {n} summary rows for {len(signals)} signals")


if __name__ == "__main__":
    main()
