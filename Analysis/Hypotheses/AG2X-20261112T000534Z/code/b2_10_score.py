#!/usr/bin/env python3
"""Stage B2 scoring (PRESPEC_AG2X_stageB2.md). Accessibility scorers only (ATAC, DNASE, SC_ATAC).

Phase 1: the 51 tier 2 (DIG1 ALT) loci, members and background, all 9 Stage B conditions.
Phase 2: tier 3 loci in a seeded random locus order (seed 20261111), bulk liver and hepatocyte only, in batches
of whole loci (about 1,000 variants); both conditions of a batch are scored before the next batch.
Variants scored in Stage B are not re-called (their Stage B scores are reused).
Cap: package total 230,000 hosted calls, 152,156 used before Stage B2. Main scoring stops before a chunk would
pass 229,400 (500 kept for resends near the end, 100 for the repeat). End: 100 phase 1 variants rescored in
bulk liver (cap 230,000). Every raw row carries candidate_not_adopted, package version, condition ID, call time.
"""
from __future__ import annotations

import json
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import ag2x_lib as L  # noqa: E402

B2 = L.PKG / "stageB2"
RAW, IN = B2 / "raw", B2 / "inputs"
CAP = 230_000
USED_BEFORE_B2 = 152_156
MAIN_CAP = CAP - 600
SCORERS = ["ATAC", "DNASE", "SC_ATAC"]
CONDITIONS = {
    "bulk_liver": "UBERON:0002107",
    "hepatocyte": "UBERON:0002107_CL:0000182",
    "monocyte_ctrl": "UBERON:0000178_CL:0001054",
    "neuron_ctrl": "UBERON:0000955_CL:0000540",
    "stellate": "UBERON:0002107_CL:0000632",
    "kupffer": "UBERON:0002107_CL:0000091",
    "liver_endothelial": "UBERON:0002107_CL:0000115",
    "cholangiocyte": "UBERON:0002107_CL:1000488",
    "hepg2": "EFO:0001187",
}
PHASE2_CONDITIONS = ["bulk_liver", "hepatocyte"]
BATCH = 1000


def calls_on_disk() -> int:
    n = 0
    for f in RAW.glob("*/chunk_log.jsonl"):
        n += sum(json.loads(x)["n_variants"] for x in open(f))
    return n


def plan():
    v = pd.read_csv(IN / "variants_b2.tsv", sep="\t")
    loci = pd.read_csv(IN / "loci_b2.tsv", sep="\t")
    t3 = sorted(loci.loc[loci.tier == 3, "locus_id"])
    order = list(np.random.default_rng(L.SEED).permutation(t3))
    rank = {lid: i for i, lid in enumerate(order)}
    v["first_locus"] = v.loci.str.split(";").map(lambda ls: min(ls, key=lambda x: rank.get(x, -1)))
    v["phase"] = np.where(v.tier == 2, 1, 2)
    new = v[~v.scored_in_stageB].copy()
    bad = L.assert_ref(new.uid.tolist())
    new = new[~new.uid.isin([b[0] for b in bad])]
    p1 = sorted(new.loc[new.phase == 1, "uid"])
    p2 = new[new.phase == 2].assign(r=lambda d: d.first_locus.map(rank)).sort_values(["r", "uid"])
    batches, cur, n = [], [], 0
    for lid, g in p2.groupby("r", sort=True):
        cur.append(g)
        n += len(g)
        if n >= BATCH:
            batches.append(pd.concat(cur)); cur, n = [], 0
    if cur:
        batches.append(pd.concat(cur))
    pd.DataFrame({"locus_id": order, "order": range(len(order))}).to_csv(IN / "tier3_locus_order.tsv", sep="\t", index=False)
    rows = [{"uid": u, "phase": 1, "batch": -1} for u in p1]
    for i, b in enumerate(batches):
        rows += [{"uid": u, "phase": 2, "batch": i, "locus_id": lid} for u, lid in zip(b.uid, b.first_locus)]
    pd.DataFrame(rows).to_csv(IN / "scoring_order_b2.tsv", sep="\t", index=False)
    budget = MAIN_CAP - USED_BEFORE_B2
    p1_calls = len(p1) * len(CONDITIONS)
    p2_full = len(p2) * len(PHASE2_CONDITIONS)
    p2_batches_fit = 0
    left = budget - p1_calls
    for b in batches:
        c = len(b) * len(PHASE2_CONDITIONS)
        if c > left:
            break
        left -= c
        p2_batches_fit += 1
    fit = pd.concat(batches[:p2_batches_fit]) if p2_batches_fit else p2.iloc[:0]
    info = {"ref_mismatch": len(bad), "phase1_variants": len(p1), "phase1_calls": p1_calls,
            "phase2_variants": int(len(p2)), "phase2_calls_if_all": int(p2_full), "phase2_batches": len(batches),
            "phase2_batches_expected_to_fit": p2_batches_fit, "phase2_variants_expected": int(len(fit)),
            "phase2_loci_expected": int(fit.first_locus.nunique()) if len(fit) else 0,
            "tier3_loci": len(t3), "budget_main": budget, "repeat_calls": 100,
            "planned_calls_total": int(p1_calls + len(fit) * 2 + 100),
            "variants_reused_from_stageB": int(v.scored_in_stageB.sum())}
    json.dump(info, open(IN / "scoring_plan_b2.json", "w"), indent=1)
    L.log(json.dumps(info))
    return p1, batches


def main():
    from alphagenome.models.v2 import variant_scorers as vs
    missing = [k for k in SCORERS if k not in vs.RECOMMENDED_VARIANT_SCORERS]
    if missing:
        raise SystemExit(f"scorer keys missing: {missing}")
    p1, batches = plan()
    if "--plan-only" in sys.argv:
        L.log("plan only: no hosted call made")
        return
    counter = {}
    used = USED_BEFORE_B2 + calls_on_disk()

    def keep(t):
        return t.assign(candidate_not_adopted=True)

    client = L.create_v2()
    for name in CONDITIONS:
        L.score_chunked(client, "v2", p1, CONDITIONS[name], SCORERS, RAW / f"p1_{name}", chunk=50, max_workers=4,
                        counter=counter, keep=keep, cap=MAIN_CAP, used_before=used)
        L.log(f"phase 1 {name} done; counter {counter}")
    stopped = None
    for i, b in enumerate(batches):
        uids = sorted(b.uid)
        need = len(uids) * len(PHASE2_CONDITIONS)
        spent = used + counter.get("attempted_calls", 0) + counter.get("resent_calls", 0)
        done = all((RAW / f"p2_b{i:03d}_{c}" / "chunk_log.jsonl").exists() and
                   sum(json.loads(x)["n_variants"] for x in open(RAW / f"p2_b{i:03d}_{c}" / "chunk_log.jsonl")) == len(uids)
                   for c in PHASE2_CONDITIONS)
        if not done and spent + need > MAIN_CAP:
            stopped = {"batch": i, "spent": spent, "need": need}
            L.log(f"cap reached before phase 2 batch {i}; stopping cleanly")
            break
        for c in PHASE2_CONDITIONS:
            L.score_chunked(client, "v2", uids, CONDITIONS[c], SCORERS, RAW / f"p2_b{i:03d}_{c}", chunk=50,
                            max_workers=4, counter=counter, keep=keep, cap=MAIN_CAP, used_before=used)
        L.log(f"phase 2 batch {i} done; counter {counter}")
    rep = sorted(np.random.default_rng(L.SEED).choice(p1, size=min(100, len(p1)), replace=False).tolist())
    L.score_chunked(client, "v2", rep, CONDITIONS["bulk_liver"], SCORERS, RAW / "repeat_bulk_liver", chunk=50,
                    max_workers=4, counter=counter, keep=keep, cap=CAP, used_before=used)
    summary = {"utc": L.utc(), "package_version": L.package_version(), "this_job_counter": counter,
               "stopped_before_batch": stopped, "hosted_calls_stageB2_on_disk": calls_on_disk(),
               "package_total_completed": USED_BEFORE_B2 + calls_on_disk(),
               "package_total_with_this_job_resends": USED_BEFORE_B2 + calls_on_disk() + counter.get("resent_calls", 0)}
    json.dump(summary, open(RAW / "scoring_summary.json", "w"), indent=1)
    L.log(json.dumps(summary))


if __name__ == "__main__":
    main()
