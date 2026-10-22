"""Independent re-derivation for the F2d blind re-check.

Own implementation (does not import the producer or H58 code). Reads release
inputs only. Bootstrap uses its own RNG stream (seed 20260923, plain
default_rng), so intervals are expected to agree with H58 within Monte Carlo
error, not bit-for-bit.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
OUT = Path(__file__).resolve().parent.parent
RUN = ROOT / "figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816"
COORDS = RUN / "multiplicity/corrected_locus_coordinates.tsv"
BASE = ROOT / "GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321"
MEMBERS = BASE / "credible_set_members_annotated.tsv"
H58 = ROOT / "Analysis/Hypotheses/H58-20261015T112423Z/jobs/5800/H58-20261015T113519Z_22252009"
SITES = H58 / "site_reannotated.tsv"
COLOC = ROOT / "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv"
GENE = ROOT / "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv"
PAIRS = ROOT / "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_signal_pairs.tsv.gz"
TIER = ROOT / "GWAS/finemapping/config/gwas_trait_tier.tsv"
SEED = 20260923
B = 10_000
DIR, ENZ = "tier1_direct_masld_pdff", "tier2_liver_enzyme"
CODING = ["protein_altering_pip_mass", "canonical_splice_pip_mass"]
ARCH = CODING + ["synonymous_or_utr_pip_mass", "other_noncoding_pip_mass"]
res = {}


def clusters(df, gap=1_000_000):
    df = df.sort_values(["chrom", "pos"]).copy()
    cid = []
    for ch, blk in df.groupby("chrom", sort=False):
        new = np.r_[True, np.diff(blk.pos.to_numpy()) > gap]
        cid.append(pd.Series([f"{ch}:{i}" for i in np.cumsum(new)], index=blk.index))
    df["cl"] = pd.concat(cid).reindex(df.index)
    return df


def cl_table(df):
    return df.groupby(["trait_scope", "cl"], as_index=False)[["num", "den"]].mean()


def frac(t, scope):
    s = t[t.trait_scope == scope]
    return s.num.sum() / s.den.sum()


# ------------------------------------------------------------- coding
a = pd.read_csv(COORDS, sep="\t", dtype={"locus": str})
a = a[a.trait_scope.isin([DIR, ENZ])]
gate = a.pip_sum_gate_passed.astype(str).str.lower() == "true"
res["gate_fail_sets"] = a[~gate].groupby("trait_scope").size().to_dict()
a = a[gate].copy()
a["chrom"] = a.chrom.astype(str)
a["num"] = a[CODING].sum(axis=1)
a["den"] = a[ARCH].sum(axis=1)
res["n_sets"] = a.groupby("trait_scope").size().to_dict()
res["sets_by_anc"] = {k: v for k, v in a.groupby(["trait_scope", "ancestry"]).size().items().__iter__()} if False else \
    {f"{k[0]}|{k[1]}": int(v) for k, v in a.groupby(["trait_scope", "ancestry"]).size().items()}

p = clusters(a)
t = cl_table(p)
pres = t.pivot_table(index="cl", columns="trait_scope", values="num", aggfunc="size").notna()
res["pooled_clusters"] = dict(direct=int(pres[DIR].sum()), enzyme=int(pres[ENZ].sum()),
                              physical=len(pres), shared=int((pres[DIR] & pres[ENZ]).sum()))
res["pooled_frac"] = dict(direct=frac(t, DIR), enzyme=frac(t, ENZ))
res["pooled_ratio"] = res["pooled_frac"]["direct"] / res["pooled_frac"]["enzyme"]
dcl = p[p.trait_scope == DIR]
res["direct_clusters_with_EUR_direct_set"] = int(dcl[dcl.ancestry == "EUR"].cl.nunique())
res["direct_clusters_total"] = int(dcl.cl.nunique())

# concentration over direct cluster means
dt = t[t.trait_scope == DIR].copy()
dt["f"] = dt.num / dt.den
dt = dt.sort_values("f", ascending=False)
res["top5_share_of_summed_direct_cluster_fraction"] = float(dt.f.head(5).sum() / dt.f.sum())
res["top1_share"] = float(dt.f.iloc[0] / dt.f.sum())
top1 = dcl[dcl.cl == dt.cl.iloc[0]]
res["top1_cluster_sets"] = top1[["study", "chrom", "pos"]].astype(str).values.tolist()
big = dcl.groupby("cl").size().sort_values(ascending=False)
res["largest_direct_cluster_n_sets"] = int(big.iloc[0])
res["largest_direct_cluster_mean"] = float(dt.set_index("cl").loc[big.index[0], "f"])
res["largest_direct_cluster_span"] = [str(dcl[dcl.cl == big.index[0]].chrom.iloc[0]),
                                      int(dcl[dcl.cl == big.index[0]].pos.min()),
                                      int(dcl[dcl.cl == big.index[0]].pos.max())]
fe = frac(t, ENZ)
drop = []
for k in range(1, 6):
    keep = dt.iloc[k:]
    drop.append(round(float(keep.num.sum() / keep.den.sum() / fe), 4))
res["ratio_after_dropping_top_k_direct_clusters"] = drop

# subclasses on pooled clustering
p["sub"] = np.where(p.trait == "PDFF", "PDFF", np.where(p.trait_scope == DIR, "NAFLD_NASH", "ENZ"))
for sub in ["PDFF", "NAFLD_NASH"]:
    q = p[(p["sub"] == sub)]
    tt = cl_table(q)
    res[f"sub_{sub}"] = dict(n_sets=len(q), n_clusters=int(q.cl.nunique()), frac=frac(tt, DIR),
                             studies=sorted(q.study.unique()))


# stratified fractions (within-stratum clustering)
def strat(df):
    tt = cl_table(clusters(df))
    return frac(tt, DIR), frac(tt, ENZ), tt


st = {}
for s in ["EUR", "AMR"]:
    fd, fz, _ = strat(a[a.ancestry == s])
    st[s] = (fd, fz)
res["EUR"] = dict(direct=st["EUR"][0], enzyme=st["EUR"][1], diff=st["EUR"][0] - st["EUR"][1])
res["AMR"] = dict(direct=st["AMR"][0], enzyme=st["AMR"][1], diff=st["AMR"][0] - st["AMR"][1])
w = {"EUR": 16 / 26, "AMR": 10 / 26}
res["EUR+AMR"] = dict(direct=sum(w[s] * st[s][0] for s in w), enzyme=sum(w[s] * st[s][1] for s in w))
res["EUR+AMR"]["diff"] = res["EUR+AMR"]["direct"] - res["EUR+AMR"]["enzyme"]
fd, fz, _ = strat(a[a.ancestry != "EUR"])
res["nonEUR"] = dict(direct=fd, enzyme=fz, diff=fd - fz)

# ---- independent bootstrap: units = pooled clusters holding EUR or AMR rows
rng = np.random.default_rng(SEED)
pid = p.set_index("credible_set_uid").cl
rows = []
for s in ["EUR", "AMR"]:
    sub = clusters(a[a.ancestry == s].copy())
    sub["pid"] = sub.credible_set_uid.map(pid)
    assert sub.groupby("cl").pid.nunique().max() == 1
    g = sub.groupby(["trait_scope", "cl"], as_index=False).agg(num=("num", "mean"), den=("den", "mean"), pid=("pid", "first"))
    g["s"] = s
    rows.append(g)
G = pd.concat(rows)
units = sorted(G.pid.unique())
res["boot_units_EUR_AMR"] = len(units)
ui = {u: i for i, u in enumerate(units)}
keys = [(s, x) for s in ["EUR", "AMR"] for x in [DIR, ENZ]]
NUM = np.zeros((len(units), 4)); DEN = np.zeros((len(units), 4))
for r in G.itertuples():
    k = keys.index((r.s, r.trait_scope)); NUM[ui[r.pid], k] += r.num; DEN[ui[r.pid], k] += r.den
idx = rng.integers(0, len(units), size=(B, len(units)))
W = np.zeros((B, len(units)))
np.add.at(W, (np.repeat(np.arange(B), len(units)), idx.ravel()), 1)
with np.errstate(invalid="ignore", divide="ignore"):
    F = (W @ NUM) / (W @ DEN)
dE = F[:, 0] - F[:, 1]; dA = F[:, 2] - F[:, 3]; dC = w["EUR"] * dE + w["AMR"] * dA
q = lambda x, lo, hi: [float(np.nanpercentile(x, lo)), float(np.nanpercentile(x, hi))]
L = 100 * 0.05 / 6
res["boot"] = {n: dict(ci95=q(d, 2.5, 97.5), ci98_333=q(d, L, 100 - L), sd=float(np.nanstd(d, ddof=1)),
                       mde=float(3.236 * np.nanstd(d, ddof=1)))
               for n, d in [("EUR", dE), ("AMR", dA), ("EUR+AMR", dC)]}


def boot_two(df, nb=B, seed_off=1):
    """physical-cluster bootstrap of direct-enzyme difference, clusters on df."""
    tt = cl_table(clusters(df))
    wide = tt.pivot(index="cl", columns="trait_scope", values=["num", "den"]).fillna(0.0)
    r = np.random.default_rng([SEED, seed_off])
    ii = r.integers(0, len(wide), size=(nb, len(wide)))
    nd, dd = wide[("num", DIR)].to_numpy()[ii].sum(1), wide[("den", DIR)].to_numpy()[ii].sum(1)
    ne, de = wide[("num", ENZ)].to_numpy()[ii].sum(1), wide[("den", ENZ)].to_numpy()[ii].sum(1)
    with np.errstate(invalid="ignore", divide="ignore"):
        d = nd / dd - ne / de
    return q(d, 2.5, 97.5), len(wide)


res["boot_pooled_all"] = boot_two(a, seed_off=1)
res["boot_nonEUR"] = boot_two(a[a.ancestry != "EUR"], seed_off=2)
# subclass bootstrap over pooled physical clusters (subset rows + all enzyme)
for j, sub in enumerate(["PDFF", "NAFLD_NASH"]):
    keep = p[(p["sub"] == sub) | (p["sub"] == "ENZ")]
    tt = keep.groupby(["trait_scope", "cl"], as_index=False)[["num", "den"]].mean()
    wide = tt.pivot(index="cl", columns="trait_scope", values=["num", "den"]).fillna(0.0)
    r = np.random.default_rng([SEED, 10 + j])
    ii = r.integers(0, len(wide), size=(B, len(wide)))
    with np.errstate(invalid="ignore", divide="ignore"):
        d = wide[("num", DIR)].to_numpy()[ii].sum(1) / wide[("den", DIR)].to_numpy()[ii].sum(1) - \
            wide[("num", ENZ)].to_numpy()[ii].sum(1) / wide[("den", ENZ)].to_numpy()[ii].sum(1)
    res[f"sub_{sub}"]["diff"] = frac(tt, DIR) - frac(tt, ENZ)
    res[f"sub_{sub}"]["ci95"] = q(d, 2.5, 97.5)

# ------------------------------------------------------------- context
mem = pd.read_csv(MEMBERS, sep="\t", low_memory=False, dtype={"locus": str})
mem = mem[mem.credible_set_uid.isin(a.credible_set_uid)].copy()
res["gated_member_rows"] = len(mem)
mem["chromosome"] = mem.chromosome.astype(str)
cons = mem.consequence.astype(str)
utr = cons.str.contains(r"(?:^|,)[35]_prime_UTR_variant(?:$|,)")
syn = cons.str.contains(r"(?:^|,)synonymous_variant(?:$|,)")
mem["nc"] = (mem.consequence_category == "other_noncoding") | ((mem.consequence_category == "synonymous_or_utr") & utr & ~syn)
unres = (mem.promoter_context_status != "resolved") | (mem.lineage_accessibility_status != "resolved")
res["adopted_unresolved_rows_by_anc"] = {s: [int(unres[mem.ancestry == s].sum()), int((mem.ancestry == s).sum())]
                                         for s in sorted(mem.ancestry.unique())}
sites = pd.read_csv(SITES, sep="\t", dtype={"chromosome": str})
res["sites"] = dict(n=len(sites), n_map_counts=sites.n_map.value_counts().to_dict())
tf = lambda x: x.astype(str).str.upper() == "TRUE"
for c in ["promoter_new", "accessible_new", "abc_new"]:
    sites[c] = tf(sites[c])
m2 = mem.merge(sites, on=["chromosome", "position"], how="left")
assert m2.n_map.notna().all()
res["unique_member_sites"] = int(mem[["chromosome", "position"]].drop_duplicates().shape[0])
# agreement on previously resolved rows
fl = lambda col: m2[col].astype(str).str.lower() == "true"
resolved = (m2.promoter_context_status == "resolved") & (m2.lineage_accessibility_status == "resolved")
agree = (fl("promoter_proximal") == m2.promoter_new) & (fl("lineage_accessible") == m2.accessible_new) & (fl("abc_enhancer_overlap") == m2.abc_new)
res["resolved_rows"] = int(resolved.sum()); res["resolved_rows_agree"] = int((agree & resolved).sum())
nc = m2[m2.nc].copy()
pip = nc.normalized_susie_pip.astype(float)
for name, flags in {"adopted": None, "repaired": True}.items():
    if flags:
        pr, ac, ab = nc.promoter_new, nc.accessible_new, nc.abc_new
        prr = acr = pd.Series(True, index=nc.index)
    else:
        pr, ac, ab = (nc.promoter_proximal.astype(str).str.lower() == "true", nc.lineage_accessible.astype(str).str.lower() == "true",
                      nc.abc_enhancer_overlap.astype(str).str.lower() == "true")
        prr, acr = nc.promoter_context_status == "resolved", nc.lineage_accessibility_status == "resolved"
    masses = pd.DataFrame({"scope": nc.trait_scope, "nc": pip, "promoter": pip * (prr & pr), "abc": pip * ab,
                           "accessible": pip * (acr & ac), "none": pip * (prr & acr & ~pr & ~ab & ~ac),
                           "unresolved": pip * (~prr | ~acr)})
    s = masses.groupby("scope").sum()
    res[f"context_summed_{name}"] = {sc: {c: round(float(s.loc[sc, c] / s.loc[sc, "nc"]), 6) for c in
                                          ["promoter", "abc", "accessible", "none", "unresolved"]} for sc in s.index}

cs = pd.read_csv(H58 / "cause_split_noncoding_rows.tsv", sep="\t")
res["cause_split_noncoding_pip"] = cs.groupby("cause").pip.sum().round(3).to_dict()
res["cause_split_noncoding_pip_EUR"] = round(float(cs[cs.ancestry == "EUR"].pip.sum()), 3)

# ------------------------------------------------------------- COLOC
tier = pd.read_csv(TIER, sep="\t")
main = tier[(tier.placement == "main") & tier.tier.isin([1, 2])]
res["main_t12_studies_in_tier_table"] = len(main)
c = pd.read_csv(COLOC, usecols=["gwas_name", "gene", "ensembl", "PP.H4.susie", "n_cs_pairs", "eligibility_source", "method"],
                low_memory=False)
c = c.merge(main[["study_name", "tier", "trait"]], left_on="gwas_name", right_on="study_name")
res["coloc_strata_present"] = int(c.gwas_name.nunique())
res["coloc_rows_by_source_tier"] = {f"{k[0]}|{k[1]}": int(v) for k, v in c.groupby(["eligibility_source", "tier"]).size().items()}
tested = c[c["PP.H4.susie"].notna()]
g = tested.groupby("tier").agg(rows=("gene", "size"), pos=("PP.H4.susie", lambda x: int((x > 0.5).sum())),
                               cs_pairs=("n_cs_pairs", "sum"))
res["tested_nonnull_susie"] = {int(k): dict(rows=int(r.rows), pos=int(r.pos), cs_pairs=float(r.cs_pairs),
                                            rate=round(r.pos / r.rows, 4)) for k, r in g.iterrows()}
tb = [[g.loc[1, "pos"], g.loc[1, "rows"] - g.loc[1, "pos"]], [g.loc[2, "pos"], g.loc[2, "rows"] - g.loc[2, "pos"]]]
orr, pv = fisher_exact(tb)
res["fisher_nominal"] = dict(OR=float(orr), p=float(pv))
res["fold_rows"] = float(g.loc[2, "rows"] / g.loc[1, "rows"]); res["fold_cs_pairs"] = float(g.loc[2, "cs_pairs"] / g.loc[1, "cs_pairs"])
pos = tested[tested["PP.H4.susie"] > 0.5]
gd, ge = set(pos[pos.tier == 1].ensembl), set(pos[pos.tier == 2].ensembl)
res["genes_pos"] = dict(direct=len(gd), enzyme=len(ge), both=len(gd & ge), direct_only=len(gd - ge), union=len(gd | ge))
res["positive_pairs_total"] = int(len(pos))
res["max_PPH4_by_source"] = c.groupby("eligibility_source")["PP.H4.susie"].max().to_dict()
gl = pd.read_csv(GENE)
pos_g = gl[gl.coloc_best_susie_pp4 > 0.5]
res["gene_level_susie_pos"] = int(len(pos_g))
res["driving_tier"] = pos_g.driving_tier.value_counts().to_dict()
res["driving_trait"] = pos_g.driving_trait.value_counts().to_dict()
sp = pd.read_csv(PAIRS, sep="\t")
sp = sp.merge(main[["study_name", "tier"]], left_on="gwas", right_on="study_name")
el = sp[sp.eligible.astype(str).str.upper() == "TRUE"]
res["eligible_cs_pairs"] = {int(k): dict(pairs=int(len(v)), hits=int((v["PP.H4"] > 0.5).sum())) for k, v in el.groupby("tier")}
res["eligible_direct_by_study"] = el[el.tier == 1].gwas.value_counts().to_dict()
res["rerun_gene_gwas_pairs"] = int(sp[["gene", "gwas"]].drop_duplicates().shape[0])

# Fig 2A support from positive pairs
pp = pd.read_csv(H58 / "FIG2A_positive_pairs.tsv", sep="\t")
res["fig2a_positive_pairs"] = len(pp)
res["fig2a_genes"] = int(pp.gene.nunique())
res["fig2a_genes_same_250kb"] = int(pp[pp.same_250000.astype(str) == "True"].gene.nunique())
nocs = {"MVP_ALT_EUR", "MVP_AST_EUR", "MVP_NAFLD_EUR"}
only = pp.groupby("gene").gwas_name.apply(lambda s: set(s) <= nocs)
res["fig2a_genes_only_MVP_EUR"] = int(only.sum())
res["fig2a_pos_pairs_vs_coloc"] = int(len(pos))

json.dump(res, open(OUT / "rederive_recheck.json", "w"), indent=1, default=str)
print(json.dumps(res, indent=1, default=str))
