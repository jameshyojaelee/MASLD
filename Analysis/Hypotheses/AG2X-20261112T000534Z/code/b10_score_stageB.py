#!/usr/bin/env python3
"""Stage B scoring (PRESPEC_AG2X_stageB.md). Accessibility scorers only; no expression scorer at GWAS loci.

Order (so a partial run is interpretable): tier 1 in all conditions, then tier 2, then tier 3.
Within a tier, conditions in the fixed order of CONDITIONS. 50-variant chunks, 4 workers, resumable.
Hosted-call cap: 150,000 for the package; 24,610 spent in Stage 0 + A; 200 reserved for the end repeat.
If tiers 1-3 in all conditions would exceed the cap, tier 3 is scored only in bulk liver, hepatocyte,
monocyte and neuron (decided before any call from the planned count).
Every raw row carries candidate_not_adopted, package version, condition ID and call time.
"""
from __future__ import annotations

import json
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import ag2x_lib as L  # noqa: E402

B = L.PKG / "stageB"
RAW = B / "raw"
SUB = L.REPO / "GWAS/finemapping/results/seqfunc/variant_substrate_hg38.tsv"
PACKAGE_CAP = 150_000
SPENT_BEFORE_B = 24_610
REPEAT_RESERVE = 200
SCORERS = ["ATAC", "DNASE", "SC_ATAC"]
CONDITIONS = {  # fixed order
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
TIER3_REDUCED = ["bulk_liver", "hepatocyte", "monocyte_ctrl", "neuron_ctrl"]


def calls_on_disk() -> int:
    n = 0
    for f in RAW.glob("*/chunk_log.jsonl"):
        n += sum(json.loads(x)["n_variants"] for x in open(f))
    return n


def fix_indels(v):
    """Substrate indels carry no hg38 allele columns; rebuild them from the hg19 alleles at pos_hg38 and keep the
    orientation whose REF matches the FASTA. Returns the table with fixed uids and the list of failures."""
    s = pd.read_csv(SUB, sep="\t", dtype=str)
    s = s[s.is_indel == "TRUE"]
    bad = v.uid.str.endswith("::")
    out, fail = {}, []
    for u in v.uid[bad]:
        c, p = u.split(":")[:2]
        r = s[(s.chr == c) & (s.pos_hg38 == p)]
        done = False
        for _, x in r.iterrows():
            for ref, alt in ((x.ref, x.alt), (x.alt, x.ref)):
                cand = f"{c}:{p}:{ref}:{alt}"
                if not L.assert_ref([cand]):
                    out[u] = cand
                    done = True
                    break
            if done:
                break
        if not done:
            fail.append(u)
    pd.DataFrame({"uid_b0": list(out), "uid": list(out.values())}).to_csv(B / "inputs/indel_uid_map.tsv", sep="\t",
                                                                          index=False)
    v = v.copy()
    v["uid"] = v.uid.map(lambda x: out.get(x, x))
    return v[~v.uid.isin(fail)], fail, len(out)


def main():
    counter = {}
    v = pd.read_csv(B / "inputs/variants_to_score.tsv", sep="\t", dtype={"tier": int})
    v[["roles", "loci"]] = v[["roles", "loci"]].fillna("")
    v, indel_fail, indel_fixed = fix_indels(v)
    # SORT1 rs12740374 is the DIG1_6 lead (tier 2); it is moved to tier 1 so it is scored first, role "control" added
    sort1 = "chr1:109274968:G:T"
    if sort1 not in set(v.uid):
        v = pd.concat([v, pd.DataFrame({"uid": [sort1], "tier": [2], "roles": [""], "loci": [""]})], ignore_index=True)
    v.loc[v.uid == sort1, "tier"] = 1
    v.loc[v.uid == sort1, "roles"] = v.loc[v.uid == sort1, "roles"] + ";control_SORT1"
    # merge rows that became the same uid after the indel rebuild: lowest tier, roles and loci joined
    v = (v.groupby("uid", as_index=False)
          .agg(tier=("tier", "min"), roles=("roles", lambda x: ";".join(sorted(set(";".join(x).split(";")) - {""}))),
               loci=("loci", lambda x: ";".join(sorted(set(";".join(x).split(";")) - {""})))))
    hsd = [u for u in v.uid if u.startswith("chr4:87310240:")]
    bad = L.assert_ref(v.uid.tolist())
    v = v[~v.uid.isin([b[0] for b in bad])]
    v.to_csv(B / "inputs/variants_scored_final.tsv", sep="\t", index=False)
    n = v.tier.value_counts().to_dict()
    full = sum(n.values()) * len(CONDITIONS) + 1  # +1 splice call for the HSD17B13 control
    budget = PACKAGE_CAP - SPENT_BEFORE_B - REPEAT_RESERVE
    tier3_conditions = list(CONDITIONS) if full <= budget else TIER3_REDUCED
    planned = (n.get(1, 0) + n.get(2, 0)) * len(CONDITIONS) + n.get(3, 0) * len(tier3_conditions) + 1
    plan = {"variants_by_tier": {int(k): int(x) for k, x in n.items()}, "indels_rebuilt": indel_fixed,
            "indels_failed": indel_fail, "ref_mismatch": len(bad), "planned_calls": planned, "budget": budget,
            "tier3_conditions": tier3_conditions, "hsd17b13_uid": hsd}
    json.dump(plan, open(B / "inputs/scoring_plan.json", "w"), indent=1)
    L.log(json.dumps(plan))
    if planned > budget:
        raise SystemExit("planned calls exceed the budget; stopping before any call")
    from alphagenome.models.v2 import variant_scorers as vs
    missing = [k for k in SCORERS + ["SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS"]
               if k not in vs.RECOMMENDED_VARIANT_SCORERS]
    if missing:
        raise SystemExit(f"scorer keys not in v2 RECOMMENDED_VARIANT_SCORERS: {missing}")
    if "--plan-only" in sys.argv:
        L.log("plan only: no hosted call made")
        return

    used = SPENT_BEFORE_B + calls_on_disk()
    cap = PACKAGE_CAP - REPEAT_RESERVE

    def keep(t):
        return t.assign(candidate_not_adopted=True)

    client = L.create_v2()
    for tier in (1, 2, 3):
        uids = sorted(v.loc[v.tier == tier, "uid"])
        conds = list(CONDITIONS) if tier < 3 else tier3_conditions
        for name in conds:
            L.score_chunked(client, "v2", uids, CONDITIONS[name], SCORERS, RAW / f"t{tier}_{name}", chunk=50,
                            max_workers=4, counter=counter, keep=keep, cap=cap, used_before=used)
            L.log(f"tier {tier} {name} done; counter {counter}")
    # HSD17B13 rs72613567 splice scorers in bulk liver (descriptive control; splice not benchmarked)
    if hsd:
        try:
            L.score_chunked(client, "v2", hsd[:1], CONDITIONS["bulk_liver"],
                            ["SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS"], RAW / "control_hsd17b13_splice",
                            chunk=1, max_workers=1, counter=counter, keep=keep, cap=cap, used_before=used)
        except Exception as e:  # descriptive control only; a failure is logged and does not stop the repeat check
            L.log(f"HSD17B13 splice control failed: {type(e).__name__}")
    # end-of-run repeat: 200 tier 1 variants, seeded, bulk liver
    rng = np.random.default_rng(L.SEED)
    t1 = sorted(v.loc[v.tier == 1, "uid"])
    rep = sorted(rng.choice(t1, size=min(200, len(t1)), replace=False).tolist())
    L.score_chunked(client, "v2", rep, CONDITIONS["bulk_liver"], SCORERS, RAW / "repeat_bulk_liver", chunk=50,
                    max_workers=4, counter=counter, keep=keep, cap=PACKAGE_CAP, used_before=used)
    summary = {"utc": L.utc(), "package_version": L.package_version(), "this_job_counter": counter,
               "hosted_calls_stageB_on_disk": calls_on_disk(),
               "hosted_calls_package_total": SPENT_BEFORE_B + calls_on_disk()}
    json.dump(summary, open(RAW / "scoring_summary.json", "w"), indent=1)
    L.log(json.dumps(summary))


if __name__ == "__main__":
    main()
