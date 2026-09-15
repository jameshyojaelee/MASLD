#!/usr/bin/env python3
"""B1 post-hoc checks over the response layer, written before the layer's numbers were read.

C1  Every signal's rank-1 variant equals the top_variant that posterior_summaries recorded for it.
C2  A_direct and B_direct share 10 of 19 analysis blocks, so the spec forbids a naive contrast between
    them; this prints the observed sharing rather than assuming it.
C3  Independent 1-Mb blocks per universe and per trait_class, so a later test knows its resampling n.
C4  The readable share by universe, coverage state and variant class.
C5  Measured channels: rows in a donor snATAC peak, rows with a DA evidence state, rows with a measured
    allelic-imbalance value, rows with an eQTL direction, and the indels the model API rescued.
C6  No unserved row carries an Atlas quantile (the layer's own guard, re-run from the written file).

Usage: b1_checks.py <response layer directory>
"""

from __future__ import annotations

import csv
import gzip
import json
import pathlib
import sys
from collections import defaultdict

PROJECT = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
TRACK0 = PROJECT / "GWAS/finemapping/results/alphagenome_atlas/run-20260909T153939Z/tables"
ATLAS_QUANTILES = ["atac_liver_quantile", "dnase_liver_quantile", "h3k27ac_liver_quantile",
                   "rna_target_gene_quantile", "rna_strongest_abs_quantile",
                   "splice_site_usage_liver_quantile", "avi_quantile"]


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

    # C2
    blocks = defaultdict(set)
    for r in rows:
        blocks[r["universe"]].add(r["analysis_block"])
    out["C2_analysis_block_sharing"] = {
        "blocks_per_universe": {u: len(b) for u, b in sorted(blocks.items())},
        "A_direct_and_B_direct_shared": len(blocks["A_direct"] & blocks["B_direct"]),
        "A_direct_only": len(blocks["A_direct"] - blocks["B_direct"]),
        "B_direct_only": len(blocks["B_direct"] - blocks["A_direct"]),
        "direct_and_C_enzyme_shared": len((blocks["A_direct"] | blocks["B_direct"]) & blocks["C_enzyme"])}

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

    json.dump(out, (d / "b1_checks.json").open("w"), indent=1, default=float)
    print(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
