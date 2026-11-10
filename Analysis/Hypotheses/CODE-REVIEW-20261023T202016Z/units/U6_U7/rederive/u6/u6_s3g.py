#!/usr/bin/env python3
"""Code review U6: Figure S3G matched pairs, re-derived from the legend and Results text (Methods give no procedure):
426 genetic genes (SuSiE PP.H4 > 0.5) in the 14,704 stable-ID joint universe, each matched 1:1 without replacement to a
non-genetic gene on bulk AveExpr (greedy nearest neighbour in a random visiting order; 0.2-SD caliper per the producer
header). Calls: TREAT interval null, treat_fdr < 0.05 (lfc 0.25). Paired OR = genetic-only / control-only;
exact McNemar = two-sided binomial on discordant pairs.
Arms: pool A = all non-genetic joint genes (as in the printed panel); pool B = SuSiE tested_le_0.5 genes only
(genes with no multi-signal COLOC result or a SuSiE-untestable pair excluded, as Figure 1B does).
Same-seed attempts use np.random.default_rng(20260817) draws 1..3 as a pandas sample random_state with two case orders.
Usage: u6_s3g.py <out_dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats

OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=True)
ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
b = pd.read_csv(ROOT / "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv")
c = pd.read_csv(ROOT / "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv")
b["id"] = b.gene.str.split(".").str[0]; c["id"] = c.ensembl.astype(str).str.split(".").str[0]
j = b.merge(c, on="id", validate="one_to_one")            # bulk file order kept
j["gen"] = j.susie_state_t12 == "positive"
j["call"] = j.treat_fdr < 0.05
prod = pd.read_csv(ROOT / "figures/candidates/fig3g-coloc-r2-r8-20261009T171017Z/source_tables/fig3g_expression_matched_pairs.tsv", sep="\t")
lines = []
def say(s): print(s); lines.append(s)
say(f"joint {len(j)} genetic {int(j.gen.sum())} TREAT calls in joint {int(j.call.sum())} genetic TREAT {int((j.gen & j.call).sum())}")

def greedy(cases, pool, order, cal):
    """nearest available pool gene within caliper, visiting cases in `order`; returns pool index per case (-1 = none)"""
    pv = pool.AveExpr.to_numpy(); srt = np.argsort(pv, kind="stable"); sv = pv[srt]
    used = np.zeros(len(pv), bool); cv = cases.AveExpr.to_numpy(); out = np.full(len(cv), -1)
    for i in order:
        k = np.searchsorted(sv, cv[i]); lo, hi = k - 1, k
        while lo >= 0 and used[srt[lo]]: lo -= 1
        while hi < len(sv) and used[srt[hi]]: hi += 1
        cand = [(abs(sv[x] - cv[i]), srt[x]) for x in (lo, hi) if 0 <= x < len(sv)]
        if not cand: continue
        d, p = min(cand)
        if d <= cal: out[i] = p; used[p] = True
    return out
def summ(cases, pool, m):
    ok = m >= 0; g = cases.call.to_numpy()[ok]; k = pool.call.to_numpy()[m[ok]]
    a, bb, cc, d = int((g & ~k).sum()), int((~g & k).sum()), int((~g & ~k).sum()), int((g & k).sum())
    p = stats.binomtest(a, a + bb, 0.5).pvalue if a + bb else np.nan
    return dict(pairs=int(ok.sum()), gen_only=a, ctrl_only=bb, neither=cc, both=d, OR=a / bb if bb else np.inf, p=p)

sd = j.AveExpr.std(ddof=1); cal = 0.2 * sd
say(f"AveExpr SD over joint {sd:.4f}; caliper 0.2 SD = {cal:.4f}")
pools = {"A_all_nongenetic": j[~j.gen], "B_tested_le_0.5": j[j.susie_state_t12 == "tested_le_0.5"]}
cases_by_order = {"bulk_order": j[j.gen].reset_index(drop=True), "id_sorted": j[j.gen].sort_values("id").reset_index(drop=True)}
rows = []
# same-seed attempts on pool A
for oname, cases in cases_by_order.items():
    rng = np.random.default_rng(20260817)
    for draw in (1, 2, 3):
        rs = int(rng.integers(1 << 31))
        order = cases.sample(frac=1.0, random_state=rs).index.to_numpy()
        pool = pools["A_all_nongenetic"].reset_index(drop=True)
        m = greedy(cases, pool, order, cal)
        s = summ(cases, pool, m)
        mine = pd.DataFrame({"g": cases.id, "c": [pool.id.iat[x] if x >= 0 else None for x in m]})
        same = mine.merge(prod, left_on=["g", "c"], right_on=["genetic_gene", "matched_control_gene"]).shape[0]
        say(f"same seed 20260817 draw {draw} case order {oname}: {s} | pairs identical to printed table {same}/426")
        rows.append(dict(arm="A", seed=20260817, draw=draw, order=oname, **s, identical=same))
# printed pair table recomputed from the adopted bulk file
pt = prod.merge(j[["id", "call", "AveExpr"]].rename(columns={"id": "genetic_gene", "call": "g_call", "AveExpr": "g_ave"}), on="genetic_gene") \
         .merge(j[["id", "call", "AveExpr"]].rename(columns={"id": "matched_control_gene", "call": "c_call", "AveExpr": "c_ave"}), on="matched_control_gene")
say(f"printed pair table vs adopted bulk file: genetic calls agree {int((pt.g_call == pt.genetic_treat_supported).sum())}/426, control calls agree {int((pt.c_call == pt.matched_control_treat_supported).sum())}/426, "
    f"max |dAveExpr| within pair {(pt.g_ave - pt.c_ave).abs().max():.4f}; counts from table: gen_only {int((pt.g_call & ~pt.c_call).sum())} ctrl_only {int((~pt.g_call & pt.c_call).sum())} both {int((pt.g_call & pt.c_call).sum())}")
# alternate seeds and 1,000-order distribution, both pools
for arm, pool in pools.items():
    pool = pool.reset_index(drop=True); cases = cases_by_order["bulk_order"]
    for seed in (20261023, 20261024):
        rng = np.random.default_rng(seed); m = greedy(cases, pool, rng.permutation(len(cases)), cal); s = summ(cases, pool, m)
        say(f"{arm} seed {seed}: {s}"); rows.append(dict(arm=arm, seed=seed, draw=0, order="perm", **s, identical=np.nan))
    rng = np.random.default_rng(20261025); res = pd.DataFrame([summ(cases, pool, greedy(cases, pool, rng.permutation(len(cases)), cal)) for _ in range(1000)])
    say(f"{arm} 1,000 visiting orders (seed 20261025): pairs {res.pairs.min()}-{res.pairs.max()}; gen_only {res.gen_only.min()}-{res.gen_only.max()} (median {res.gen_only.median()}); "
        f"ctrl_only {res.ctrl_only.min()}-{res.ctrl_only.max()} (median {res.ctrl_only.median()}); OR {res.OR.min():.3f}-{res.OR.max():.3f} (median {res.OR.median():.3f}); "
        f"McNemar p {res.p.min():.4f}-{res.p.max():.4f}; share p<0.05 {(res.p < 0.05).mean():.3f}")
    res.to_csv(OUT / f"s3g_1000orders_{arm}.tsv", sep="\t", index=False)
pd.DataFrame(rows).to_csv(OUT / "s3g_runs.tsv", sep="\t", index=False)
(OUT / "s3g_log.txt").write_text("\n".join(lines) + "\n")
