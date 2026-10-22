# F2d repair re-derivation (2026-10-15). Point values only; no new bootstrap.
# A: ancestry composition of the 6,415 gated credible sets and 233 pooled clusters
#    (checker failure R02-01 "16 of 19").
# B: COLOC r2 denominators per tested gene-trait row and per tested credible-set pair
#    (checker failure R14-03).
# C: noncoding-context estimators after the H58 position repair (text problem on
#    fig2.md:107 "weighting each credible set equally").
# Inputs are the adopted release files; the H58 re-annotated site table is the only
# H58 product used (C), and the producer's clustering is imported read-only.
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
H58 = ROOT / "Analysis/Hypotheses/H58-20261015T112423Z/jobs/5800/H58-20261015T113519Z_22252009"
sys.path.insert(0, str(H58 / "scripts"))
import h58_common as h  # noqa: E402  read-only import of loaders
from h58_common import D, E, m  # noqa: E402

OUT = Path(sys.argv[1])
res = {}

# ---------------------------------------------------------------- A
a = h.load_coding()
assert len(a) == 6415
a = m.assign_clusters(a.copy())
res["n_sets"] = a.trait_scope.value_counts().to_dict()
cl_scope = a.groupby("cluster_id").trait_scope.agg(lambda s: set(s))
res["clusters_total"] = int(len(cl_scope))
res["clusters_direct"] = int(sum(D in s for s in cl_scope))
res["clusters_enzyme"] = int(sum(E in s for s in cl_scope))
res["clusters_shared"] = int(sum(len(s) == 2 for s in cl_scope))
dd = a[a.trait_scope == D]
ee = a[a.trait_scope == E]
res["direct_sets_by_ancestry"] = dd.ancestry.value_counts().to_dict()
res["enzyme_sets_by_ancestry"] = ee.ancestry.value_counts().to_dict()
dcl = set(dd.cluster_id)
res["direct_clusters_with_EUR_direct_set"] = int(dd[dd.ancestry == "EUR"].cluster_id.nunique())
res["direct_clusters_with_any_EUR_set"] = int(a[(a.cluster_id.isin(dcl)) & (a.ancestry == "EUR")].cluster_id.nunique())
res["direct_clusters_with_only_EUR_direct_sets"] = int(
    dd.groupby("cluster_id").ancestry.agg(lambda s: set(s) == {"EUR"}).sum())
res["direct_clusters_with_nonEUR_direct_set"] = int(dd[dd.ancestry != "EUR"].cluster_id.nunique())
eur_d = dd[dd.ancestry == "EUR"].copy()
eur_d = m.assign_clusters(eur_d)
res["EUR_direct_clusters_within_EUR_only"] = int(eur_d.cluster_id.nunique())
# pooled fractions (producer estimator: per-cluster means, ratio of sums)
t = m.cluster_table(a, "coding_mass", "arch_total")
fr = t.groupby("trait_scope").apply(lambda g: g.num.sum() / g.den.sum())
res["coding_fraction_pooled"] = {k: round(float(v), 6) for k, v in fr.items()}
# direct subclasses against all enzyme clusters (pooled clustering)
for lab, mask in [("PDFF", dd.trait == "PDFF"), ("NAFLD_NASH", dd.trait != "PDFF")]:
    sub = dd[mask]
    tt = m.cluster_table(sub, "coding_mass", "arch_total")
    res[f"subclass_{lab}"] = dict(n_sets=int(len(sub)), n_clusters=int(sub.cluster_id.nunique()),
                                  fraction=round(float(tt.num.sum() / tt.den.sum()), 6),
                                  studies=sorted(sub.study.unique().tolist()))

# ---------------------------------------------------------------- B
tier = pd.read_csv(ROOT / "GWAS/finemapping/config/gwas_trait_tier.tsv", sep="\t")
main = tier[(tier.placement == "main") & tier.tier.isin([1, 2])].set_index("study_name")
assert len(main) == 35
c = pd.read_csv(h.COLOC, usecols=["gwas_name", "gene", "ensembl", "PP.H4.susie", "n_cs_pairs",
                                  "n_cs_pairs_eligible", "eligibility_source", "method"],
                low_memory=False)
c = c[c.gwas_name.isin(main.index)].copy()
c["tier"] = c.gwas_name.map(main.tier)
res["rows_by_source_tier"] = {f"{k[0]}|{k[1]}": int(v) for k, v in
                              c.groupby(["eligibility_source", "tier"]).size().items()}
tested = c[c["PP.H4.susie"].notna()].copy()
tested["pos"] = tested["PP.H4.susie"] > 0.5
g = tested.groupby("tier").agg(rows=("pos", "size"), pos=("pos", "sum"),
                               cs_pairs=("n_cs_pairs", "sum"))
res["tested_rows"] = {int(k): dict(rows=int(r.rows), pos=int(r.pos), cs_pairs=float(r.cs_pairs),
                                   rate=round(r.pos / r.rows, 4)) for k, r in g.iterrows()}
f = fisher_exact([[g.loc[1, "pos"], g.loc[1, "rows"] - g.loc[1, "pos"]],
                  [g.loc[2, "pos"], g.loc[2, "rows"] - g.loc[2, "pos"]]])
res["tested_rows_fisher_nominal"] = dict(OR=round(float(f[0]), 4), p=float(f[1]))
nr = tested[tested.eligibility_source == "not_rerun_negative"]
res["max_PPH4_not_rerun"] = float(nr["PP.H4.susie"].max())
rr = c[c.eligibility_source == "rerun_shared_posterior_check"]
res["rerun_rows_by_tier"] = {int(k): dict(rows=int(len(v)), pos=int((v["PP.H4.susie"] > 0.5).sum()))
                             for k, v in rr.groupby("tier")}
res["fold_tested_rows_enzyme_over_direct"] = round(g.loc[2, "rows"] / g.loc[1, "rows"], 2)
res["fold_tested_cs_pairs_enzyme_over_direct"] = round(g.loc[2, "cs_pairs"] / g.loc[1, "cs_pairs"], 2)
# genes
gpos = tested[tested.pos]
g1, g2 = set(gpos[gpos.tier == 1].ensembl), set(gpos[gpos.tier == 2].ensembl)
res["genes_pos"] = dict(direct=len(g1), enzyme=len(g2), both=len(g1 & g2),
                        direct_only=len(g1 - g2), union=len(g1 | g2))
gt1, gt2 = set(tested[tested.tier == 1].ensembl), set(tested[tested.tier == 2].ensembl)
res["genes_tested"] = dict(direct=len(gt1), enzyme=len(gt2), union=len(gt1 | gt2))
# signal pairs
sp = pd.read_csv(h.COLOC.parent / "susie_coloc_signal_pairs.tsv.gz", sep="\t")
sp = sp[sp.gwas.isin(main.index)].copy()
sp["tier"] = sp.gwas.map(main.tier)
el = sp.eligible.astype(str).str.upper().isin(["TRUE", "1"])
sp = sp[el]
sp["hit"] = sp["PP.H4"] > 0.5
res["eligible_pairs"] = {int(k): dict(pairs=int(len(v)), hits=int(v.hit.sum()))
                         for k, v in sp.groupby("tier")}
res["eligible_direct_pairs_by_study"] = sp[sp.tier == 1].gwas.value_counts().to_dict()
gl = pd.read_csv(h.COLOC.parent / "gene_level_coloc_tier12.csv")
hh = gl[gl.coloc_best_susie_pp4 > 0.5]
res["driving_tier"] = hh.driving_tier.value_counts().to_dict()

# ---------------------------------------------------------------- C
mem = h.load_gated_members()
sites = h.load_sites(H58)
rep = h.repaired_members(mem, sites)
ctx_rep = h.context_per_set(rep)
ctx_adp = h.context_per_set(mem)
X = ["promoter_mass", "abc_mass", "accessible_mass", "resolved_no_context_mass",
     "unresolved_context_mass"]
out_c = {}
for lab, frame in [("adopted", ctx_adp), ("repaired", ctx_rep)]:
    s = frame.groupby("trait_scope")[X + ["noncoding_pip_mass"]].sum()
    ratio_of_sums = {sc: {x: round(float(s.loc[sc, x] / s.loc[sc, "noncoding_pip_mass"]), 6)
                          for x in X} for sc in (D, E)}
    fr_ = frame[frame.noncoding_pip_mass > 0].copy()
    eq = {sc: {x: round(float((fr_[fr_.trait_scope == sc][x] /
                               fr_[fr_.trait_scope == sc].noncoding_pip_mass).mean()), 6)
               for x in X} for sc in (D, E)}
    cw = {}
    cfr = m.assign_clusters(frame.copy())
    for x in X:
        tt = m.cluster_table(cfr, x, "noncoding_pip_mass")
        cw[x] = {sc: round(float(tt[tt.trait_scope == sc].num.sum() /
                                 tt[tt.trait_scope == sc].den.sum()), 6) for sc in (D, E)}
    out_c[lab] = dict(summed_over_sets=ratio_of_sums, mean_of_set_fractions=eq,
                      cluster_weighted=cw,
                      n_sets_with_noncoding={sc: int((fr_.trait_scope == sc).sum()) for sc in (D, E)})
res["context"] = out_c
# unresolved rows by ancestry before repair (noncoding member rows)
nc = mem.copy()
res["adopted_unresolved_rows_by_ancestry"] = {}
if "promoter_context_status" in nc.columns:
    for anc, v in nc.groupby("ancestry"):
        res["adopted_unresolved_rows_by_ancestry"][anc] = dict(
            rows=int(len(v)), unresolved=int((v.promoter_context_status != "resolved").sum()))

(OUT / "rederive.json").write_text(json.dumps(res, indent=1, default=str))
print(json.dumps(res, indent=1, default=str))
