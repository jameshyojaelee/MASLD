"""Repair-round Figure 7E table for review round 4, workstream F7a.

Starts again from the release table (fig7-source-exposure-correction-20261014/tables), not
from the first-round copy.
  R08-04  7E treatment response: diet-arm MDE as threshold, statistic oriented so that
          positive is the prespecified direction. Rationale reworded: the point estimate
          is opposite to the prespecified direction and its interval spans zero; the
          source analysis records the result as negative and underpowered.
  R08-06  7E signed eQTL direction: the released row is kept as released (tested
          negative). The prespecified rule 2 tests the label-reliability ceiling, and every
          stratum passes it. The k=1 sequence-attributable bound is reported as a
          sensitivity in derivation_checks.tsv only.
  R08-05  7F two-layer range re-derived from the release gain table (no panel change).
"""
import csv
import json
import math
import os

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
EX = f"{ROOT}/Analysis/MASLD_Model_Benchmark/executions"
REL = f"{EX}/fig7-source-exposure-correction-20261014/tables"
RD = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = f"{RD}/exec/tables"
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
diff = wad["diff"]
mde_diet = wad["mde_80pct"]
lo, hi = wad["ci95"]
mde_re = Z80 * wad["sd_of_deltas"] * math.sqrt(1 / wad["n_resolved"] + 1 / wad["n_persisted"])
oriented = -diff
checks += [
    ("R08-04 diet arm n resolved / persisted", f"{wad['n_resolved']} / {wad['n_persisted']}"),
    ("R08-04 diet diff resolved minus persisted", diff),
    ("R08-04 diet 95% CI (resolved minus persisted)", f"{lo:.6f} to {hi:.6f}"),
    ("R08-04 diet permutation p", wad["permutation_p"]),
    ("R08-04 diet MDE80 file", mde_diet),
    ("R08-04 diet MDE80 recomputed 2.802*sd*sqrt(1/6+1/16)", mde_re),
    ("R08-04 oriented observed (persisted minus resolved)", oriented),
    ("R08-04 ratio oriented / diet MDE", oriented / mde_diet),
    ("R08-04 oriented 95% CI / diet MDE", f"{-hi / mde_diet:.4f} to {-lo / mde_diet:.4f}"),
    ("R08-04 CI excludes a correct-direction effect of one diet MDE", lo > -mde_diet),
]
assert abs(mde_re - mde_diet) < 1e-3
assert lo < 0 < hi  # interval spans zero: the point estimate alone is not a reversed effect

f7e, v = read_tsv(f"{REL}/panel_7e_verdicts.tsv")
hit = 0
for r in v:
    if r["question"] == "severity scores move correctly when the same person is treated":
        hit += 1
        assert r["verdict"] == "tested_negative"  # released verdict, unchanged
        r["statistic"] = ("within-participant delta, persisted minus resolved, diet arm only "
                          "(positive = the prespecified direction)")
        r["observed"] = repr(oriented)
        r["threshold"] = repr(mde_diet)
        r["threshold_name"] = "80% MDE of the diet-arm contrast (6 resolved vs 16 persisted)"
        r["decision_rule"] = (
            "prespecified direction (PRESPECIFICATION.md, written for the pooled P2 contrast): "
            "the score falls further where NASH resolves, and a reversed direction is reported as "
            "a negative result. The diet arm is the only arm with more than one persister (surgery "
            "14 vs 1). There, resolved minus persisted = +0.0593 (95% CI -0.0261 to +0.1456; "
            "permutation p 0.232): the point estimate is opposite to the prespecified direction and "
            "the interval spans zero, so it is not a reversed effect either. The interval excludes "
            "a correct-direction effect as large as the diet-arm MDE (-0.138). The source analysis "
            "records the result as negative and underpowered (RESULTS.md). The earlier threshold "
            "0.1001 was the MDE of the pooled 20 vs 17 contrast.")
        r["source"] = ("Analysis/MASLD_Model_Benchmark/executions/gse83452-intervention-response-"
                       "20260831/results_20260831/within_arm_decisive.json")
    if r["question"] == "signed cis-eQTL direction from DNA-LM embeddings":
        hit += 1
        assert r["verdict"] == "tested_negative"  # kept as released (R08-06 repair)
assert hit == 2
write_tsv(f"{TAB}/panel_7e_verdicts.tsv", f7e, v)

# ------------------------------------------------------------------ R08-06 sensitivity only
strata = json.load(open(f"{EX}/eqtl_sign_reliability-20260831T195702Z/logs/03_strata_intervals.json"))
checks.append(("R08-06 03_strata_intervals.json top-level keys", ",".join(list(strata)[:12])))
f7d, b = read_tsv(f"{REL}/panel_7d_bounds.tsv")
BA = [r for r in b if r["block"].startswith("a.")]
for s in dict.fromkeys(r["row"] for r in BA):
    g = {r["quantity"]: float(r["value"]) for r in BA if r["row"] == s}
    half = (g["achieved AUC 95% CI high"] - g["achieved AUC 95% CI low"]) / 2
    mde = Z80 * half / 1.96
    bound = g["sequence-attributable bound"] - g["measured chance"]
    checks.append((f"R08-06 {s}: k=1 bound-chance / CI half-width / MDE80 / bound<half / bound<MDE",
                   f"{bound:.4f} / {half:.4f} / {mde:.4f} / {bound < half} / {bound < mde}"))

# ------------------------------------------------------------------ R08-05 range, release table
_, g7f = read_tsv(f"{REL}/panel_7f_gain.tsv")
two = [r for r in g7f if r["block"].startswith("b.") and r["head_id"] == "two_layer_mlp"]
for r in two:
    r["model_id"] = r["model_id"].split("::")[0]  # model_id carries a "::head" suffix
pre = [r for r in two if r["eligibility"].startswith("interpretable") and r["model_id"] != "scvi_liver_latent"
       or r["model_id"] == "transcriptformer_tf_sapiens"]
for r in two:
    checks.append((f"R08-05 two-layer {r['model_id']}: gain [lo, hi] / clean / eligibility",
                   f"{float(r['gain']):+.5f} [{float(r['lo']):+.5f}, {float(r['hi']):+.5f}] / "
                   f"{r['clean_studies']} / {r['eligibility'][:60]}"))
pg = sorted(float(r["gain"]) for r in pre)
checks.append(("R08-05 two-layer range, pretrained encoders with >=1 clean study (excl. scVI liver, Geneformer)",
               f"{pg[0]:+.5f} to {pg[-1]:+.5f} over {len(pg)} encoders"))
allg = sorted(float(r["gain"]) for r in two if r["model_id"] != "geneformer_v2_316m")
checks.append(("R08-05 two-layer range incl. scVI liver (excl. Geneformer)", f"{allg[0]:+.5f} to {allg[-1]:+.5f}"))

with open(f"{RD}/exec/derivation_checks.tsv", "w") as f:
    f.write("check\tvalue\n")
    for k, val in checks:
        f.write(f"{k}\t{val}\n")
for k, val in checks:
    print(k, "|", val)
