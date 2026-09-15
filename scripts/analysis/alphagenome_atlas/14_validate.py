#!/usr/bin/env python3
"""Step 14: fail-closed reconciliation of the produced tables (protocol §10 checks)."""

from __future__ import annotations

import csv
import json
import math
import pathlib
from collections import defaultdict

import pandas as pd

import lib_atlas as la

P = la.prespec()
ROOT = la.out_root()
TABLES = ROOT / "tables"
RAW = ROOT / "raw"
OUT = ROOT / "validation"


def refuse_if_written(out: pathlib.Path) -> None:
    """Write-once guard. An empty directory made by job setup is fine; one that already holds a result is not."""
    if out.exists() and any(out.iterdir()):
        raise la.ContractError(f"refusing to overwrite: {out}")
    out.mkdir(parents=True, exist_ok=True)


def check(cond: bool, msg: str, results: list) -> None:
    results.append({"check": msg, "pass": bool(cond)})
    if not cond:
        la.log(f"FAIL: {msg}")


def main() -> None:
    refuse_if_written(OUT)
    res = []
    sig = la.read_tsv(TABLES / "eligible_signals.tsv")
    n = defaultdict(int)
    for s in sig:
        n[s["universe"]] += 1
    check(n["A_direct"] == 29 and n["C_enzyme"] == 787, f"universe counts A=29 C=787 (got {dict(n)})", res)
    check(len({s["signal_uid"] for s in sig}) == len(sig), "signal_uid unique", res)
    check(all(s["posterior_definition"] for s in sig), "posterior_definition present on every signal", res)
    # weights: one row per (signal, variant); mass ≤ 1 + tol for A/C
    mass, seen, dup = defaultdict(float), set(), 0
    with la.open_text(TABLES / "signal_variant_weights.tsv.gz") as h:
        for r in csv.DictReader(h, delimiter="\t"):
            k = (r["signal_uid"], r["source_variant_id"])
            dup += k in seen; seen.add(k)
            if r["mapping_status"] == "mapped":
                mass[r["signal_uid"]] += float(r["weight"])
    check(dup == 0, "no duplicated (signal, variant) weight rows", res)
    uni = {s["signal_uid"]: s["universe"] for s in sig}
    check(all(m <= 1 + 1e-6 for k, m in mass.items() if uni[k] != "B_direct"), "COLOC mapped mass ≤ 1 per signal", res)
    # archive count = query set
    qs = sum(1 for _ in open(TABLES / "query_set.tsv")) - 1
    archived = set()
    for stage in ("atlas_direct", "atlas_enzyme"):
        for c in (RAW / stage).glob("chunk_*"):
            archived.update(json.load((c / "request.json").open())["variants"])
    check(len(archived) == qs, f"archived unique variants ({len(archived)}) == query set ({qs})", res)
    avail = pd.read_csv(TABLES / "atlas_availability.tsv.gz", sep="\t")
    check(len(avail) == qs, f"availability rows ({len(avail)}) == query set", res)
    # posterior summaries: C_liver ≤ C_atlas ≤ C_query ≤ C_map ≤ 1+tol ; signed NaN when unsigned
    ps = pd.read_csv(TABLES / "posterior_summaries.tsv.gz", sep="\t")
    tol = 1e-6
    check(bool(((ps["C_liver"] <= ps["C_atlas"] + tol) & (ps["C_atlas"] <= ps["C_query"] + tol) & (ps["C_query"] <= ps["C_map"] + tol)).all()), "coverage layers nested", res)
    unsigned = ps[~ps["is_signed"]]
    check(bool(unsigned["S_liver_median"].isna().all()), "unsigned scorers carry no signed summary", res)
    # hand re-derivation of one signal's C and S for RNA_SEQ target gene
    liver = pd.read_csv(TABLES / "liver_summaries.tsv.gz", sep="\t", keep_default_na=False)
    target_rows = ps[(ps["scorer"] == "RNA_SEQ") & (ps["gene_scope"] == "target_gene") & (ps["C_liver"] > 0)].sort_values("signal_uid")
    if len(target_rows):
        r = target_rows.iloc[0]
        w = {}
        with la.open_text(TABLES / "signal_variant_weights.tsv.gz") as h:
            for x in csv.DictReader(h, delimiter="\t"):
                if x["signal_uid"] == r["signal_uid"] and x["in_query_set"] == "True" and x["mapping_status"] == "mapped":
                    w[x["variant_uid"]] = float(x["weight"])
        g = liver[(liver["scorer"] == "RNA_SEQ") & (liver["gene_id"] == r["gene_id"]) & (liver["variant_uid"].isin(w))]
        sc = dict(zip(g["variant_uid"], g["liver_median_raw"].astype(float)))
        c_hand = sum(w[v] for v in sc)
        s_hand = sum(w[v] * sc[v] for v in sc)
        check(abs(c_hand - r["C_liver"]) < 1e-9 and abs(s_hand - r["S_liver_median"]) < 1e-9, f"hand re-derivation of C and S for {r['signal_uid']} (C {c_hand:.6f} vs {r['C_liver']:.6f})", res)
    # profiles cover every signal
    prof = pd.read_csv(TABLES / "signal_profiles.tsv", sep="\t")
    check(len(prof) == len(sig), "signal_profiles has one row per signal", res)
    cen = pd.read_csv(TABLES / "chromatin_census.tsv", sep="\t")
    check(len(cen) == len(sig), "chromatin_census has one row per signal", res)
    # scorer metadata used, signed flags from server
    md = la.read_tsv(TABLES / "scorer_metadata.tsv")
    check(len(md) == 22, f"22 scorers in metadata (got {len(md)})", res)
    gate = json.load((TABLES / "ag_atlas_gate_replication.json").open()) if (TABLES / "ag_atlas_gate_replication.json").exists() else {}
    check("auroc" in gate, "gate replication computed", res)
    pd.DataFrame(res).to_csv(OUT / "checks.tsv", sep="\t", index=False)
    ok = all(r["pass"] for r in res)
    la.log(f"VALIDATION {'PASS' if ok else 'FAIL'}: {sum(r['pass'] for r in res)}/{len(res)} checks")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
