#!/usr/bin/env python3
"""B1 v2 post-hoc checks over the response layer.

C1  Every signal's rank-1 variant equals the top_variant that posterior_summaries recorded for it.
C2  Analysis-block sharing between universes, with 'no analysis block' counted as MISSING, not as a block,
    and with the 1-Mb block overlaps printed separately because they are a different unit.
C3  Independent 1-Mb blocks per universe and per trait_class, so a later test knows its resampling n.
C4  The readable share by universe, coverage state and variant class.
C5  Measured channels, each with its composition: rows in a donor snATAC peak, rows with a DA evidence
    state (and which states), rows with a measured allelic-imbalance value, rows with an eQTL direction,
    and the indels the model API rescued (with the 1-Mb block recovered from the hg38 identity).
C6  No unserved row carries an Atlas quantile (the layer's own guard, re-run from the written file).
C7  Explanation classes per signal.
C8  Every signed column names the allele or the group its sign refers to, on every row that carries a value.
C9  The weight-sum guard, stated on the full posterior vector and on the exported (floored) rows separately.
C10 Allelic-imbalance donor counts per SITE with the site n, beside the per-row figure that is not the unit.
C11 Distinct hg38 uids in the Track 0 crosswalk, and how many carry both allele orders.

Usage: b1_checks_v2.py <response layer directory>
"""

from __future__ import annotations

import csv
import gzip
import json
import pathlib
import statistics
import sys
from collections import defaultdict

PROJECT = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
TRACK0 = PROJECT / "GWAS/finemapping/results/alphagenome_atlas/run-20260909T153939Z/tables"
ATLAS_QUANTILES = ["atac_primary_liver_quantile", "dnase_primary_liver_quantile",
                   "h3k27ac_primary_liver_quantile", "rna_target_gene_primary_liver_quantile",
                   "rna_strongest_abs_primary_liver_quantile",
                   "splice_site_usage_primary_liver_quantile", "avi_quantile"]
# every signed value column -> the column that states what its sign refers to
SIGNED_COLUMNS = {
    "ase_gse281367_mean_log2_alt_over_ref": "ase_sign_refers_to",
    "ase_gse244832_mean_log2_alt_over_ref": "ase_sign_refers_to",
    "gwas_beta_allele1": "direction_sign_refers_to",
    "eqtl_beta_allele1": "direction_sign_refers_to",
    "pred_allele1_liver_rna": "direction_sign_refers_to",
    "model_api_rna_log2": "model_api_sign_refers_to",
    "model_api_splice_log2": "model_api_sign_refers_to",
    "model_api_atac_log2": "model_api_sign_refers_to",
    "model_api_dnase_log2": "model_api_sign_refers_to",
    "model_api_h3k27ac_log2": "model_api_sign_refers_to",
    "measured_atac_logFC_gse244832": "measured_atac_logFC_gse244832_sign_refers_to",
    "measured_atac_logFC_gse281367": "measured_atac_logFC_gse281367_sign_refers_to",
}


def rows_of(path):
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as h:
        yield from csv.DictReader(h, delimiter="\t")


def main():
    d = pathlib.Path(sys.argv[1])
    rows = list(rows_of(d / "response_layer_variants.tsv.gz"))
    sigs = list(rows_of(d / "response_layer_signals.tsv"))
    out = {}

    # C1
    top_ps = {}
    for r in rows_of(TRACK0 / "posterior_summaries.tsv.gz"):
        top_ps.setdefault(r["signal_uid"], set()).add(r["top_variant"])
    mine = {r["signal_uid"]: r["variant_uid"] for r in rows if r["weight_rank"] == "1"}
    agree = sum(1 for s, v in mine.items() if s in top_ps and v in top_ps[s])
    absent = [s for s in mine if s not in top_ps]
    out["C1_rank1_matches_posterior_summaries_top_variant"] = {
        "n_signals_checked": len(mine) - len(absent), "n_agree": agree,
        "n_signals_absent_from_posterior_summaries": len(absent),
        "n_signals_with_more_than_one_top_variant_across_scorers":
            sum(1 for s in top_ps if len(top_ps[s]) > 1 and s in mine),
        "disagreements": [{"signal_uid": s, "b1_rank1": v, "posterior_summaries": sorted(top_ps[s])}
                          for s, v in mine.items() if s in top_ps and v not in top_ps[s]][:20]}

    # C2. 'no analysis block' is missing, not a block: counting the empty string adds a phantom block to
    # every universe that has one, and makes two universes look like they share a locus when they do not.
    blocks, mb = defaultdict(set), defaultdict(set)
    nb_signals = defaultdict(set)
    for r in rows:
        if r["analysis_block"]:
            blocks[r["universe"]].add(r["analysis_block"])
        else:
            nb_signals[r["universe"]].add(r["signal_uid"])
        if r["block_1mb"]:
            mb[r["universe"]].add(r["block_1mb"])
    out["C2_analysis_block_sharing"] = {
        "unit": "Track 0 analysis_block; the 1-Mb block is a different unit and is reported beside it",
        "named_blocks_per_universe": {u: len(b) for u, b in sorted(blocks.items())},
        "signals_with_NO_analysis_block_per_universe": {u: len(s) for u, s in sorted(nb_signals.items())},
        "blocks_per_universe_if_empty_counted_as_a_block":
            {u: len(b) + (1 if nb_signals.get(u) else 0) for u, b in sorted(blocks.items())},
        "A_direct_and_B_direct_shared_ANALYSIS_blocks": len(blocks["A_direct"] & blocks["B_direct"]),
        "A_direct_only": len(blocks["A_direct"] - blocks["B_direct"]),
        "B_direct_only": len(blocks["B_direct"] - blocks["A_direct"]),
        "direct_and_C_enzyme_shared_ANALYSIS_blocks":
            len((blocks["A_direct"] | blocks["B_direct"]) & blocks["C_enzyme"]),
        "A_direct_and_B_direct_shared_1MB_blocks": len(mb["A_direct"] & mb["B_direct"]),
        "direct_and_C_enzyme_shared_1MB_blocks":
            len((mb["A_direct"] | mb["B_direct"]) & mb["C_enzyme"])}

    # C3
    def blockcount(sel, key):
        g = defaultdict(set)
        for r in sel:
            if r["block_1mb"]:
                g[r[key]].add(r["block_1mb"])
        return {k: len(v) for k, v in sorted(g.items())}
    out["C3_independent_1mb_blocks"] = {
        "total": len({r["block_1mb"] for r in rows if r["block_1mb"]}),
        "by_universe": blockcount(rows, "universe"),
        "by_trait_class": blockcount(rows, "trait_class"),
        "by_posterior_definition": blockcount(rows, "posterior_definition"),
        "readable_rows_only_by_universe": blockcount([r for r in rows if r["row_readable"] == "True"], "universe")}

    # C4
    def frac(sel):
        return {"n": len(sel), "n_readable": sum(r["row_readable"] == "True" for r in sel),
                "share": (sum(r["row_readable"] == "True" for r in sel) / len(sel)) if sel else None}
    out["C4_readable_share"] = {
        "all": frac(rows),
        "by_universe": {u: frac([r for r in rows if r["universe"] == u])
                        for u in sorted({r["universe"] for r in rows})},
        "by_coverage_state": {s: frac([r for r in rows if r["signal_atlas_coverage_state"] == s])
                              for s in sorted({r["signal_atlas_coverage_state"] for r in rows})},
        "by_variant_class": {c: frac([r for r in rows if r["variant_class"] == c])
                             for c in sorted({r["variant_class"] for r in rows})}}

    # C5
    def blocks_of(sel):
        return len({r["block_1mb"] for r in sel if r["block_1mb"]})
    ch = {
        "in_donor_snATAC_peak": [r for r in rows if r["measured_peak"]],
        "da_state_supported_or_source_dependent": [
            r for r in rows if set(r["measured_da_evidence_state"].split(";")) & {"supported", "source_dependent"}],
        "measured_allelic_imbalance": [r for r in rows if r["ase_gse281367_mean_log2_alt_over_ref"]
                                       or r["ase_gse244832_mean_log2_alt_over_ref"]],
        "eqtl_direction": [r for r in rows if r["eqtl_beta_allele1"]],
        "model_api_indel_rescue": [r for r in rows if r["model_api_splice_log2"]],
        "atlas_served": [r for r in rows if r["atlas_served"] == "True"]}
    out["C5_measured_channels"] = {k: {"n_rows": len(v), "n_signals": len({r["signal_uid"] for r in v}),
                                       "n_1mb_blocks": blocks_of(v),
                                       "n_signals_by_universe": {u: len({r["signal_uid"] for r in v
                                                                         if r["universe"] == u})
                                                                 for u in ("A_direct", "B_direct", "C_enzyme")}}
                                   for k, v in ch.items()}
    da = ch["da_state_supported_or_source_dependent"]
    out["C5_measured_channels"]["da_state_supported_or_source_dependent"]["composition"] = {
        s: sum(1 for r in da if r["measured_da_evidence_state"] == s)
        for s in sorted({r["measured_da_evidence_state"] for r in da})}
    out["C5_measured_channels"]["da_state_supported_or_source_dependent"]["n_rows_exactly_supported"] = sum(
        1 for r in rows if r["measured_da_evidence_state"] == "supported")
    api = ch["model_api_indel_rescue"]
    out["C5_measured_channels"]["model_api_indel_rescue"]["blocks_recovered"] = {
        "n_rows_with_empty_block_1mb": sum(1 for r in api if not r["block_1mb"]),
        "n_rows_with_a_recovered_block": sum(1 for r in api if r["block_1mb_recovered"]),
        "n_distinct_recovered_1mb_blocks": len({r["block_1mb_recovered"] for r in api
                                                if r["block_1mb_recovered"]}),
        "recovered_from": {k: sum(1 for r in api if r["block_1mb_recovered_from"] == k)
                           for k in sorted({r["block_1mb_recovered_from"] for r in api})}}

    # C6
    bad = [(r["signal_uid"], r["variant_uid"], c) for r in rows if r["atlas_served"] == "False"
           for c in ATLAS_QUANTILES if str(r.get(c, "")).strip()]
    out["C6_no_unserved_row_carries_an_atlas_quantile"] = {"n_violations": len(bad), "examples": bad[:5]}

    # explanation classes
    out["C7_signal_classes"] = {
        c: {v: sum(s[c] == v for s in sigs) for v in sorted({s[c] for s in sigs})}
        for c in ("class_coding", "class_local_regulatory", "class_rna_processing", "class_candidate_long_range")}
    out["C7_signal_classes"]["locus_state"] = {
        v: sum(s["locus_state"] == v for s in sigs) for v in sorted({s["locus_state"] for s in sigs})}
    out["C7_signal_classes"]["by_universe"] = {
        u: {v: sum(1 for s in sigs if s["universe"] == u and s["locus_state"] == v)
            for v in sorted({s["locus_state"] for s in sigs})}
        for u in ("A_direct", "B_direct", "C_enzyme")}

    # C8: every signed column states what its sign refers to, on every row that carries a value.
    cols = set(rows[0].keys())
    c8 = {}
    for col, ocol in SIGNED_COLUMNS.items():
        if col not in cols:
            c8[col] = {"present": False}
            continue
        vals = [r for r in rows if str(r[col]).strip()]
        c8[col] = {"present": True, "n_values": len(vals),
                   "orientation_column": (ocol if ocol in cols else "ABSENT FROM TABLE"),
                   "n_rows_missing_the_orientation": (
                       sum(1 for r in vals if not str(r.get(ocol, "")).strip()) if ocol in cols else len(vals)),
                   "distinct_orientation_texts": sorted({r.get(ocol, "") for r in vals})[:3]}
    out["C8_signed_columns_state_their_orientation"] = {
        "per_column": c8,
        "columns_with_any_row_missing_an_orientation": sorted(
            c for c, v in c8.items() if v.get("n_rows_missing_the_orientation", 0) > 0 or not v["present"])}

    # C9: the weight guard, on the full posterior vector and on the exported rows separately.
    gate = list(rows_of(d / "coverage_gate.tsv"))
    total_mass = {g["signal_uid"]: float(g["total_mass"]) for g in gate}
    exported = defaultdict(float)
    for r in rows:
        exported[r["signal_uid"]] += float(r["weight"])
    full = defaultdict(float)
    for r in rows_of(TRACK0 / "signal_variant_weights.tsv.gz"):
        if r["signal_uid"] in total_mass:
            full[r["signal_uid"]] += float(r["weight"])
    out["C9_weight_sums"] = {
        "what_the_producer_guards": ("the FULL posterior weight vector sums to the coverage table's "
                                     "total_mass to within 1e-6, for every signal"),
        "n_signals": len(total_mass),
        "n_signals_full_vector_within_1e-6": sum(1 for s, t in total_mass.items() if abs(full[s] - t) < 1e-6),
        "max_abs_difference_full_vector": max(abs(full[s] - t) for s, t in total_mass.items()),
        "n_signals_EXPORTED_rows_within_1e-6": sum(1 for s, t in total_mass.items()
                                                   if abs(exported[s] - t) < 1e-6),
        "why_the_exported_sum_differs": ("exported rows are the variants at or above 0.01 posterior plus the "
                                         "rank-1 variant; the floor removes mass at every signal whose "
                                         "credible set has a tail")}

    # C10: allelic-imbalance donors per SITE.
    ase_rows = [r for r in rows if r["ase_gse281367_mean_log2_alt_over_ref"]
                or r["ase_gse244832_mean_log2_alt_over_ref"]]
    seen, per_site = set(), []
    for r in ase_rows:
        site = r["variant_uid"] or r["variant_key"]
        if site not in seen:
            seen.add(site)
            per_site.append(r)
    c10 = {"n_rows": len(ase_rows), "n_distinct_sites": len(per_site),
           "unit": "donor; the site is the unit of the median, not the (signal, variant) row"}
    for ds in ("gse281367", "gse244832"):
        col = f"ase_{ds}_n_het_donors"
        sv = sorted(float(r[col]) for r in per_site if r[col])
        rv = sorted(float(r[col]) for r in ase_rows if r[col])
        c10[ds] = {"n_sites": len(sv), "median_per_site": (statistics.median(sv) if sv else None),
                   "min": (sv[0] if sv else None), "max": (sv[-1] if sv else None),
                   "n_rows": len(rv), "median_per_row_NOT_the_unit": (statistics.median(rv) if rv else None)}
    out["C10_allelic_imbalance_donors_per_site"] = c10

    # C11: the crosswalk uid count that DEFECTS.md B1-D3 quotes.
    uid_all, uid_swapmix = set(), defaultdict(set)
    for r in rows_of(TRACK0 / "variant_crosswalk.tsv.gz"):
        u = r["variant_uid"]
        uid_all.add(u)
        uid_swapmix[u].add(r["allele_swap"])
    out["C11_crosswalk_uids"] = {
        "n_distinct_variant_uid_including_empty": len(uid_all),
        "n_distinct_variant_uid_non_empty": len(uid_all - {""}),
        "n_uids_carrying_both_none_and_swapped": sum(1 for u, s in uid_swapmix.items()
                                                     if u and {"none", "swapped"} <= s)}

    json.dump(out, (d / "b1_checks.json").open("w"), indent=1, default=float)
    print(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
