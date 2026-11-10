#!/usr/bin/env python3
"""Code review U7 (H71): independent T1 for the Kumak NAS and fibrosis lists, from the Methods text:
U_eval = 3,769 genes (SuSiE positive or tested_le_0.5 in the stable-ID join); genetic = SuSiE PP.H4 > 0.5;
10 equal-count AveExpr bins (ties by Ensembl ID); exact stratified hypergeometric two-sided p (|x - E| >= |O - E|);
MH OR with RBG CI; Holm over the two lists. List membership (Ensembl IDs) is taken from partA_requalification.tsv.
Usage: u7_h71.py <out_dir>"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats

OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=True)
ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
b = pd.read_csv(ROOT / "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv")
c = pd.read_csv(ROOT / "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv")
b["id"] = b.gene.str.split(".").str[0]; c["id"] = c.ensembl.astype(str).str.split(".").str[0]
j = b.merge(c, on="id", validate="one_to_one")
u = j[j.susie_state_t12.isin(["positive", "tested_le_0.5"])].copy()
u["gen"] = u.coloc_best_susie_pp4 > 0.5
u = u.sort_values(["AveExpr", "id"]).reset_index(drop=True)
u["bin"] = pd.qcut(np.arange(len(u)), 10, labels=False)
pa = pd.read_csv(ROOT / "Analysis/Hypotheses/H71-20261015T124239Z/run-20261015T125815Z/partA_requalification.tsv", sep="\t")
L = []
def say(s): print(s); L.append(s)

def t1(flag):
    d = u.assign(f=flag)
    tab = d.groupby("bin").apply(lambda g: pd.Series({"N": len(g), "K": int(g.gen.sum()), "n": int(g.f.sum()), "a": int((g.gen & g.f).sum())}))
    E = float((tab.K * tab.n / tab.N).sum()); O = int(tab.a.sum())
    pmf = np.array([1.0])
    for _, t in tab.iterrows():
        x = np.arange(0, min(t.K, t.n) + 1); pmf = np.convolve(pmf, stats.hypergeom.pmf(x, t.N, t.K, t.n))
    xs = np.arange(len(pmf)); p = float(pmf[np.abs(xs - E) >= abs(O - E) - 1e-9].sum())
    A = tab.a; B = tab.K - tab.a; C = tab.n - tab.a; D = tab.N - tab.K - tab.n + tab.a; N = tab.N
    R = A * D / N; S = B * C / N; P = (A + D) / N; Q = (B + C) / N
    OR = R.sum() / S.sum()
    v = (P * R).sum() / (2 * R.sum() ** 2) + ((P * S) + (Q * R)).sum() / (2 * R.sum() * S.sum()) + (Q * S).sum() / (2 * S.sum() ** 2)
    lo, hi = np.exp(np.log(OR) + np.array([-1, 1]) * 1.959963984540054 * np.sqrt(v))
    return dict(n_list_in_Ueval=int(flag.sum()), O=O, E=E, OE=O / E, OR=OR, lo=lo, hi=hi, p=p)

res = {}
for lst in ("Kumak_NAS", "Kumak_Fibrosis"):
    ids = set(pa[(pa.list == lst) & pa.ensembl_id.notna()].ensembl_id.str.split(".").str[0])
    inj = j[j.id.isin(ids)]
    say(f"{lst}: matched IDs {len(ids)}; in bulk and COLOC join {len(inj)}; no_susie_result {int((inj.susie_state_t12 == 'no_susie_result').sum())}")
    res[lst] = t1(u.id.isin(ids).to_numpy())
    r = res[lst]
    say(f"{lst}: in U_eval {r['n_list_in_Ueval']}; O {r['O']} E {r['E']:.3f} O/E {r['OE']:.3f}; MH OR {r['OR']:.3f} ({r['lo']:.3f} to {r['hi']:.3f}); exact p {r['p']:.4f}")
ps = np.array([res[k]["p"] for k in res]); order = np.argsort(ps); holm = np.empty(2)
holm[order] = np.maximum.accumulate(np.minimum(1, ps[order] * np.array([2, 1])))
say("Holm over the two lists: " + ", ".join(f"{k} {h:.4f}" for k, h in zip(res, holm)))
(OUT / "u7_h71_log.txt").write_text("\n".join(L) + "\n")
