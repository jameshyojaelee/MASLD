"""Build corrected Figure 7 panel tables for review round 4, workstream F7a.

Reads the producing files, never the earlier panel numbers, for every value it changes.
  R08-04  7E treatment response: diet-arm MDE as threshold, statistic oriented to the
          prespecified direction (positive = score falls more where NASH resolves).
  R08-06  7E signed eQTL direction: verdict untestable, because the sequence-attributable
          bound sits below the resolution of the test.
  R08-12  7D block b: the three transfer cohorts (achieved Spearman, 95% CI, tie ceiling).
Writes tables into ../exec/tables (copies; the originals are untouched) and a check file.
"""
import csv
import json
import math
import os

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
EX = f"{ROOT}/Analysis/MASLD_Model_Benchmark/executions"
WS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = f"{WS}/exec/tables"
Z80 = 2.802  # 1.96 + 0.842, the repository's 80%-power convention


def read_tsv(p):
    with open(p, newline="") as f:
        r = csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE)
        return r.fieldnames, list(r)


def write_tsv(p, fields, rows):
    with open(p, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", quoting=csv.QUOTE_NONE,
                           escapechar="\\", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


checks = []

# ------------------------------------------------------------------ 7E, R08-04
wad_p = f"{EX}/gse83452-intervention-response-20260831/results_20260831/within_arm_decisive.json"
wad = json.load(open(wad_p))["readouts"]["gene_median_elastic_net"]["diet_arm_only"]
diff = wad["diff"]                        # resolved minus persisted, +0.059285
mde_diet = wad["mde_80pct"]               # 0.137772
lo, hi = wad["ci95"]
# Independent recomputation of the MDE from the saved SD and group sizes.
mde_re = Z80 * wad["sd_of_deltas"] * math.sqrt(1 / wad["n_resolved"] + 1 / wad["n_persisted"])
checks.append(("R08-04 diet MDE80 file", mde_diet))
checks.append(("R08-04 diet MDE80 recomputed 2.802*sd*sqrt(1/6+1/16)", mde_re))
oriented = -diff                          # persisted minus resolved; positive = prespecified direction
checks.append(("R08-04 oriented observed", oriented))
checks.append(("R08-04 ratio oriented/diet MDE", oriented / mde_diet))
checks.append(("R08-04 old ratio diff/pooled MDE", diff / 0.10006576070998692))

f7e, v = read_tsv(f"{TAB}/panel_7e_verdicts.tsv")
n_hit = 0
for r in v:
    if r["question"] == "severity scores move correctly when the same person is treated":
        n_hit += 1
        r["statistic"] = ("within-participant delta, persisted minus resolved, diet arm only "
                          "(positive = the prespecified direction)")
        r["observed"] = repr(oriented)
        r["threshold"] = repr(mde_diet)
        r["threshold_name"] = "80% MDE of the diet-arm contrast (6 resolved vs 16 persisted)"
        r["decision_rule"] = (
            "prespecified direction: the score falls further where NASH resolves "
            "(PRESPECIFICATION.md P2 negative; a reversed direction is reported as a negative "
            "result). Diet arm resolved minus persisted = +0.0593 (95% CI -0.0261 to +0.1456; "
            "permutation p 0.232), the reversed sign. The interval excludes a correct-direction "
            "effect as large as the diet-arm MDE (-0.138), so: tested negative, direction "
            "reversed. The earlier threshold 0.1001 was the MDE of the pooled 20 vs 17 contrast.")
        r["source"] = ("Analysis/MASLD_Model_Benchmark/executions/gse83452-intervention-response-"
                       "20260831/results_20260831/within_arm_decisive.json")
    if r["question"] == "signed cis-eQTL direction from DNA-LM embeddings":
        n_hit += 1
        r["verdict"] = "untestable"

assert n_hit == 2

# ------------------------------------------------------------------ 7E, R08-06
f7d, b = read_tsv(f"{TAB}/panel_7d_bounds.tsv")
BA = [r for r in b if r["block"].startswith("a.")]
strata = []
for s in dict.fromkeys(r["row"] for r in BA):
    g = {r["quantity"]: float(r["value"]) for r in BA if r["row"] == s}
    half = (g["achieved AUC 95% CI high"] - g["achieved AUC 95% CI low"]) / 2
    mde = Z80 * half / 1.96
    bound = g["sequence-attributable bound"] - g["measured chance"]
    obs = g["achieved AUC (best embedding head)"] - g["measured chance"]
    strata.append((s, obs, bound, half, mde))
    checks.append((f"R08-06 {s}: achieved-chance / bound-chance / CI half-width / MDE80",
                   f"{obs:+.4f} / {bound:.4f} / {half:.4f} / {mde:.4f} / bound<MDE {bound < mde}"))
assert all(bd < m for _, _, bd, _, m in strata)
s0 = strata[0]
for r in v:
    if r["question"] == "signed cis-eQTL direction from DNA-LM embeddings":
        r["decision_rule"] = (
            "the largest gain any sequence model could reach (sequence-attributable bound, at "
            f"most 1.89% of held-out rows causal) is {s0[2]:.4f} above chance on all rows, below "
            f"the test's resolution half-width {s0[3]:.4f}; in every |z| stratum the bound is below "
            "the 80% MDE. A model reaching the bound would not be detected, so: untestable on "
            "these data, NOT negative. The sign-reliability ceilings (up to 1.0000) are not "
            "attainable by a sequence model and do not set the verdict.")
write_tsv(f"{TAB}/panel_7e_verdicts.tsv", f7e, v)

# ------------------------------------------------------------------ 7D, R08-12
tab_p = f"{EX}/tabicl-nonlinear-comparator-20260905T000000Z/results/results.json"
r0 = json.load(open(tab_p))["arms"]["R0_released_frozen"]["per_cohort"]
COH = [("GSE268273", "GSE268273", "fibrosis 0-4"),
       ("GSE213621_primary_3level", "GSE213621", "3-level stage"),
       ("GSE276114_primary_3level", "GSE276114", "3-level stage")]
SRC = "Analysis/MASLD_Model_Benchmark/executions/tabicl-nonlinear-comparator-20260905T000000Z/results/results.json (arms.R0_released_frozen.per_cohort)"
BLK = "b. released severity instrument in its three transfer cohorts: Spearman beside the tie-structure ceiling"
newb = []
ratios = []
for key, gse, lvl in COH:
    d = r0[key]
    row = f"{gse} {lvl}, n={d['n']} participant records"
    note = f"never in the instrument's training; ratio to ceiling {d['spearman'] / d['tie_ceiling']:.4f}"
    for q, val, kind in [("achieved Spearman (released v1.1/v1.2 weights, no refit)", d["spearman"], "achieved"),
                         ("tie-structure ceiling", d["tie_ceiling"], "ceiling"),
                         ("achieved Spearman 95% CI low", d["ci95"][0], "ci_lo"),
                         ("achieved Spearman 95% CI high", d["ci95"][1], "ci_hi")]:
        newb.append({"block": BLK, "row": row, "quantity": q, "value": repr(val), "kind": kind,
                     "note": note, "source": SRC})
    ratios.append(d["spearman"] / d["tie_ceiling"])
    checks.append((f"R08-12 {gse} n / Spearman / ceiling / ratio",
                   f"{d['n']} / {d['spearman']:.4f} / {d['tie_ceiling']:.4f} / {d['spearman'] / d['tie_ceiling']:.4f}"))
checks.append(("R08-12 ratio range", f"{min(ratios):.4f} to {max(ratios):.4f}"))
oldb = [r for r in b if r["block"].startswith("b.")]
b2 = [r for r in b if not r["block"].startswith("b.")]
idx = next(i for i, r in enumerate(b) if r["block"].startswith("b."))
b2 = b2[:idx] + newb + b2[idx:]
write_tsv(f"{TAB}/panel_7d_bounds.tsv", f7d, b2)
write_tsv(f"{WS}/exec/panel_7d_block_b_retired_rows.tsv", f7d, oldb)

with open(f"{WS}/exec/derivation_checks.tsv", "w") as f:
    f.write("check\tvalue\n")
    for k, val in checks:
        f.write(f"{k}\t{val}\n")
for k, val in checks:
    print(k, "|", val)
