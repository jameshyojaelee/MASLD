#!/usr/bin/env python3
"""Step 15 (Package C): join eligible-signal genes to the existing LeafCutter joint_v2 event family.

No new discovery run and no re-FDR among selected loci: the family-adjusted values come from the
existing joint_v2 analysis (its own analysis_contract.json). Atlas splice predictions for the signal's
variants are placed beside each measured event; a variant inside the event's intron boundaries is a
positional link only, not a validation of the variant's molecular effect (no linked genotypes).

Outputs (tables/): rna_processing_results.tsv, rna_processing_signal_summary.tsv
"""

from __future__ import annotations

import csv
import json
from collections import defaultdict

import pandas as pd

import lib_atlas as la

TABLES = la.out_root() / "tables"
J2 = la.PROJECT / "GWAS/finemapping/results/seqfunc/disease_splicing/joint_v2"
WINDOW = 0  # exact containment within intron boundaries only


def main() -> None:
    signals = [s for s in la.read_tsv(TABLES / "eligible_signals.tsv") if s["universe"] != "B_direct"]
    ev = pd.read_csv(J2 / "primary_events_annotated.tsv.gz", sep="\t", keep_default_na=False, low_memory=False)
    contract = json.load((J2 / "analysis_contract.json").open()) if (J2 / "analysis_contract.json").exists() else {}
    ev["gene_key"] = ev["gene_id"].astype(str).str.split(".").str[0]
    by_gene = {g: d for g, d in ev.groupby("gene_key")}
    liver = pd.read_csv(TABLES / "liver_summaries.tsv.gz", sep="\t", keep_default_na=False)
    junc = liver[liver["scorer"] == "SPLICE_JUNCTIONS"]
    usage = liver[liver["scorer"] == "SPLICE_SITE_USAGE"].set_index(["variant_uid", "gene_id"])
    weights = defaultdict(dict)
    with la.open_text(TABLES / "signal_variant_weights.tsv.gz") as h:
        for r in csv.DictReader(h, delimiter="\t"):
            if r["mapping_status"] == "mapped" and r["in_query_set"] == "True":
                weights[r["signal_uid"]][r["variant_uid"]] = float(r["weight"])
    rows, summary = [], []
    for s in signals:
        sig, ens = s["signal_uid"], s["ensembl"]
        d = by_gene.get(ens)
        w = weights[sig]
        if d is None or len(d) == 0:
            summary.append({"signal_uid": sig, "gene": s["gene"], "ensembl": ens, "n_events": 0, "n_events_q05": 0, "n_replicated": 0,
                            "state": "untestable_no_eligible_event_for_gene", "family": "leafcutter joint_v2 (existing)", "n_queried_variants": len(w)})
            continue
        n_q, n_rep = 0, 0
        for _, e in d.iterrows():
            try:
                b1, b2 = int(e["boundary1"]), int(e["boundary2"])
            except ValueError:
                b1 = b2 = None
            inside = [v for v in w if b1 is not None and v.split(":")[0] == str(e["chrom"]) and b1 - WINDOW <= int(v.split(":")[1]) <= b2 + WINDOW]
            inside_mass = sum(w[v] for v in inside)
            # Atlas junction matching this intron on the MEASURED convention, not a tolerance: Atlas
            # junctions are half-open and LeafCutter introns closed, so (start, end + 1) is exact.
            # Measured here: start offset 0 (17,159 hits, next best 11), end offset +1 (19,893, next best 133).
            jm = junc[(junc["gene_id"] == ens) & (junc["variant_uid"].isin(w))]
            jm = jm[[la.leafcutter_intron_key(a, b) == (b1, b2) for a, b in zip(jm["junction_start"], jm["junction_end"])]] if b1 is not None and len(jm) else jm.iloc[0:0]
            top_j = jm.reindex(jm["liver_median_raw"].abs().sort_values(ascending=False).index).head(1)
            usage_vals = [float(usage.at[(v, ens), "liver_median_quantile"]) for v in w if (v, ens) in usage.index]
            q = _f(e["cluster_q"])
            n_q += int(q == q and q < 0.05)
            n_rep += int(str(e["multicohort_replicated"]) == "True")
            rows.append({"signal_uid": sig, "gene": s["gene"], "ensembl": ens, "event_id": e["event_id"], "chrom": e["chrom"], "boundary1": e["boundary1"],
                         "boundary2": e["boundary2"], "cluster": e["cluster"], "deltapsi_Disease": e["deltapsi_Disease"], "logef_Disease": e["logef_Disease"],
                         "cluster_p": e["cluster_p"], "cluster_q_family": e["cluster_q"], "cohort_tested": e["cohort_tested"], "cohort_same_sign": e["cohort_same_sign"],
                         "loco_tested": e["loco_tested"], "multicohort_replicated": e["multicohort_replicated"], "annotated_exact": e["annotated_exact"],
                         "signal_mass_inside_intron": inside_mass, "n_signal_variants_inside_intron": len(inside),
                         "atlas_matching_junction_variant": top_j["variant_uid"].iloc[0] if len(top_j) else "",
                         "atlas_matching_junction_liver_median": top_j["liver_median_raw"].iloc[0] if len(top_j) else "",
                         "atlas_matching_junction_liver_quantile": top_j["liver_median_quantile"].iloc[0] if len(top_j) else "",
                         "atlas_splice_site_usage_max_abs_quantile": max((abs(x) for x in usage_vals), default=""),
                         "link_type": "positional_only_no_linked_genotypes", "unit": "participants (bulk cohorts); LeafCutter intron excision"})
        summary.append({"signal_uid": sig, "gene": s["gene"], "ensembl": ens, "n_events": len(d), "n_events_q05": n_q, "n_replicated": n_rep,
                        "state": "measured_events_available", "family": f"leafcutter joint_v2 (existing; contract: {contract.get('family', contract.get('analysis', 'see analysis_contract.json'))})",
                        "n_queried_variants": len(w)})
    la.write_tsv_once(TABLES / "rna_processing_results.tsv", rows, list(rows[0].keys()) if rows else ["signal_uid"])
    la.write_tsv_once(TABLES / "rna_processing_signal_summary.tsv", summary, list(summary[0].keys()))
    la.log(f"Package C: {len(summary)} signals, {len(rows)} event rows")


def _i(v):
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return -10**9


def _f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("nan")


if __name__ == "__main__":
    main()
