#!/usr/bin/env python
"""U5 independent re-derivation of Figure 2D (pooled) and its ancestry rows (H58 family C).

Written from METHODS.md (paragraph "Figure 2D, its ancestry rows, Figure S2S ...") and the input
files. Producer code was read only to learn the RNG stream (seed + crc32(fit id)) and the unit
order, which are needed for same-seed reproduction.

Inputs
  corrected_locus_coordinates.tsv  6,415 gated credible sets with restored lead positions and masses
  credible_set_pip_architecture.tsv  6,467 raw instances (gate recount)
  GWAS/finemapping/data/lead_snps/<study>_leadSNPs_tiered.tsv  lead index p
  GWAS/finemapping/results/credible_sets.csv  credible-set loci (coverage check)
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
import zlib
from pathlib import Path
from statistics import NormalDist

import numpy as np
import pandas as pd

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
COORDS = ROOT / "figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/multiplicity/corrected_locus_coordinates.tsv"
ARCH = ROOT / "GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321/credible_set_pip_architecture.tsv"
LEADS = ROOT / "GWAS/finemapping/data/lead_snps"
CREDSETS = ROOT / "GWAS/finemapping/results/credible_sets.csv"
OUT = Path(sys.argv[1])
OUT.mkdir(parents=True, exist_ok=True)

D, E = "tier1_direct_masld_pdff", "tier2_liver_enzyme"
GWS = 5e-8
GAP = 1_000_000
B = 10_000
NSPLIT = int(os.environ.get("U5_NSPLIT", "200"))  # 0 = skip calibration (login-node light check)
LO98, HI98 = 100 * 0.05 / 6, 100 - 100 * 0.05 / 6
MDE_FACTOR = NormalDist().inv_cdf(1 - 0.05 / 6) + NormalDist().inv_cdf(0.80)
LOG = []


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.append(s)


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for ch in iter(lambda: f.read(1 << 20), b""):
            h.update(ch)
    return h.hexdigest()


def gen(seed, fit_id):
    return np.random.default_rng([seed, zlib.crc32(fit_id.encode())])


# ------------------------------------------------------------------ 1. credible sets and lead filter
raw = pd.read_csv(ARCH, sep="\t", dtype={"locus": str})
raw = raw[raw.trait_scope.isin([D, E])]
gate = raw.pip_sum_gate_passed.astype(str).str.lower() == "true"
gate_num = raw.raw_cs_pip_sum.between(0.90, 1.10)
log("raw instances", len(raw), "gate flag pass", int(gate.sum()), "0.90<=sum<=1.10", int(gate_num.sum()),
    "flag==numeric", bool((gate == gate_num).all()))

cs = pd.read_csv(COORDS, sep="\t", dtype={"locus": str})
log("corrected coords rows", len(cs), "uids equal gated raw:",
    set(cs.credible_set_uid) == set(raw.loc[gate, "credible_set_uid"]))
cs["coding"] = cs.protein_altering_pip_mass + cs.canonical_splice_pip_mass
cs["den4"] = cs[["protein_altering_pip_mass", "canonical_splice_pip_mass", "synonymous_or_utr_pip_mass",
                 "other_noncoding_pip_mass"]].sum(axis=1)
log("unresolved mass >0 sets:", int((cs.unresolved_pip_mass > 0).sum()), "den4 range", cs.den4.min(), cs.den4.max())

lead_rows = []
for st in sorted(cs.study.unique()):
    f = LEADS / f"{st}_leadSNPs_tiered.tsv"
    t = pd.read_csv(f, sep="\t")
    t["study"] = st
    lead_rows.append(t)
L = pd.concat(lead_rows, ignore_index=True)
L["key"] = list(zip(L.study, L.CHR.astype(str), L.BP.astype(int)))
assert L.key.is_unique
pmap = dict(zip(L.key, L.index_p))
tmap = dict(zip(L.key, L.tier))
cs["key"] = list(zip(cs.study, cs.chrom.astype(str), cs.pos.astype(int)))
cs["lead_p"] = cs.key.map(pmap)
cs["lead_tier"] = cs.key.map(tmap)
cs["matched"] = cs.lead_p.notna()
um = ~cs.matched
log("unmatched sets", int(um.sum()), "by ancestry", cs[um].ancestry.value_counts().to_dict(),
    "by class", cs[um].trait_scope.value_counts().to_dict())
keep = (cs.matched & (cs.lead_p < GWS)) | um
k = cs[keep].copy()
counts = dict(
    n_reliable=len(cs), n_kept=len(k), n_direct=int((k.trait_scope == D).sum()), n_enzyme=int((k.trait_scope == E).sum()),
    n_unmatched_kept=int(um.sum()),
    afr_enzyme_suggestive_p_lt_gws=int(((k.ancestry == "AFR") & (k.trait_scope == E) & (k.lead_tier == "suggestive")).sum()),
    afr_suggestive_p_range=[float(k.loc[(k.lead_tier == "suggestive"), "lead_p"].min()),
                            float(k.loc[(k.lead_tier == "suggestive"), "lead_p"].max())],
    eur_direct=int(((k.ancestry == "EUR") & (k.trait_scope == D)).sum()),
    nonEUR_enzyme=int(((k.ancestry != "EUR") & (k.trait_scope == E)).sum()),
    by_class_ancestry=k.groupby(["trait_scope", "ancestry"]).size().rename("n").reset_index().to_dict("records"))
log("COUNTS", json.dumps(counts, default=int))


# ------------------------------------------------------------------ 2. clusters
def single_linkage(fr: pd.DataFrame) -> pd.Series:
    """cluster label 'chrom:k' (k 1-based within chromosome) for single linkage within GAP on lead pos."""
    lab = pd.Series(index=fr.index, dtype=object)
    for ch, g in fr.sort_values(["chrom", "pos"]).groupby("chrom", sort=True):
        newc = np.concatenate([[True], np.diff(g.pos.to_numpy()) > GAP])
        lab.loc[g.index] = [f"{ch}:{i}" for i in np.cumsum(newc)]
    return lab


k["pc"] = single_linkage(k)
ncl = k.pc.nunique()
dcl = set(k.loc[k.trait_scope == D, "pc"])
ecl = set(k.loc[k.trait_scope == E, "pc"])
eur_d_cl = set(k.loc[(k.trait_scope == D) & (k.ancestry == "EUR"), "pc"])
log("CLUSTERS total", ncl, "direct", len(dcl), "enzyme", len(ecl), "shared", len(dcl & ecl),
    "direct clusters with EUR direct set", len(dcl & eur_d_cl))
span = k.groupby("pc").pos.agg(lambda s: s.max() - s.min())
log("max cluster span bp", int(span.max()))


# ------------------------------------------------------------------ 3. pooled estimate
def class_cluster_means(fr, col="coding"):
    return fr.groupby(["trait_scope", "pc"]).agg(num=(col, "mean"), den=("den4", "mean")).reset_index()


cm = class_cluster_means(k)
frac = {x: cm.loc[cm.trait_scope == x, "num"].sum() / cm.loc[cm.trait_scope == x, "den"].sum() for x in (D, E)}
# independent plain form: mean over clusters of per-cluster mean coding
frac_plain = {x: k[k.trait_scope == x].groupby("pc").coding.mean().mean() for x in (D, E)}
diff = frac[D] - frac[E]
log(f"POOLED direct {100*frac[D]:.4f}% enzyme {100*frac[E]:.4f}% diff {diff:.6f}; plain-mean form "
    f"{100*frac_plain[D]:.4f} / {100*frac_plain[E]:.4f}")

# four-class composition (cluster-weighted)
comp = {}
for x in (D, E):
    sub = k[k.trait_scope == x]
    per = sub.groupby("pc")[["protein_altering_pip_mass", "canonical_splice_pip_mass", "synonymous_or_utr_pip_mass",
                             "other_noncoding_pip_mass"]].mean()
    comp[x] = (100 * per.sum() / per.to_numpy().sum()).round(4).to_dict()
log("COMPOSITION", json.dumps(comp))

# concentration: share of the summed direct cluster fractions carried by each direct cluster
dfrac = (cm[cm.trait_scope == D].set_index("pc").eval("num / den")).sort_values(ascending=False)
share = dfrac / dfrac.sum()
top = k[(k.trait_scope == D) & (k.pc == share.index[0])]
conc = dict(n_direct_clusters=int(len(share)), top1_share=float(share.iloc[0]), top5_share=float(share.iloc[:5].sum()),
            top1_cluster=share.index[0], top1_n_sets=int(len(top)), top1_pos=[int(p) for p in top.pos])
log("CONCENTRATION", json.dumps(conc))

# physical-cluster bootstrap. Units = pooled physical clusters; one draw takes whole clusters.
wide = cm.pivot(index="pc", columns="trait_scope", values=["num", "den"]).fillna(0.0)
units_prod = sorted(wide.index)  # producer order: lexicographic 'chrom:k'
units_mine = sorted(wide.index, key=lambda s: (int(s.split(":")[0]), int(s.split(":")[1])))


def pooled_boot(order, g):
    w = wide.loc[order]
    na, da = w[("num", D)].to_numpy(), w[("den", D)].to_numpy()
    nb, db = w[("num", E)].to_numpy(), w[("den", E)].to_numpy()
    idx = g.integers(0, len(w), size=(B, len(w)))
    with np.errstate(divide="ignore", invalid="ignore"):
        fa = na[idx].sum(1) / da[idx].sum(1)
        fb = nb[idx].sum(1) / db[idx].sum(1)
    dd = fa - fb
    ok = np.isfinite(dd)
    return dict(lo95=float(np.percentile(dd[ok], 2.5)), hi95=float(np.percentile(dd[ok], 97.5)), n_finite=int(ok.sum()),
                sd=float(np.std(dd[ok], ddof=1)))


pool_res = {
    "same_seed_20260817_producer_order": pooled_boot(units_prod, gen(20260817, "coding_mass")),
    "alt_seed_1_producer_order": pooled_boot(units_prod, gen(1, "coding_mass")),
    "alt_seed_2_producer_order": pooled_boot(units_prod, gen(2, "coding_mass")),
    "independent_order_plain_seed_20260817": pooled_boot(units_mine, np.random.default_rng(20260817)),
}
for kk, v in pool_res.items():
    log("POOLED BOOT", kk, json.dumps(v))


# ------------------------------------------------------------------ 4. ancestry family C (EUR, AMR, fixed-weight EUR+AMR)
class Strat:
    """Within-stratum clusters nested in pooled clusters; units = pooled clusters with a row in `strata`."""

    def __init__(self, fr, strata, scol, pooled_lab):
        self.keys = [(s, x) for s in strata for x in (D, E)]
        rows = []
        for s in strata:
            sub = fr[fr[scol] == s].copy()
            sub["wc"] = single_linkage(sub)
            sub["pooled"] = pooled_lab.loc[sub.index]
            assert sub.groupby("wc").pooled.nunique().max() == 1
            g = sub.groupby(["trait_scope", "wc"]).agg(num=("coding", "mean"), den=("den4", "mean"),
                                                      pooled=("pooled", "first")).reset_index()
            g["s"] = s
            rows.append(g)
        self.tab = pd.concat(rows, ignore_index=True)
        self.units = sorted(self.tab.pooled.unique())
        ui = {u: i for i, u in enumerate(self.units)}
        self.NUM = np.zeros((len(self.units), len(self.keys)))
        self.DEN = np.zeros_like(self.NUM)
        for r in self.tab.itertuples(index=False):
            j = self.keys.index((r.s, r.trait_scope))
            self.NUM[ui[r.pooled], j] += r.num
            self.DEN[ui[r.pooled], j] += r.den

    def ncl(self, s, x):
        t = self.tab
        return int(((t.s == s) & (t.trait_scope == x)).sum())

    def draws(self, g, nb=B):
        U = len(self.units)
        idx = g.integers(0, U, size=(nb, U))
        W = np.zeros((nb, U))
        np.add.at(W, (np.repeat(np.arange(nb), U), idx.ravel()), 1.0)
        den = W @ self.DEN
        with np.errstate(divide="ignore", invalid="ignore"):
            f = np.where(den > 0, (W @ self.NUM) / den, np.nan)
        return f, den

    def obs(self):
        return self.NUM.sum(0) / self.DEN.sum(0)

    def j(self, s, x):
        return self.keys.index((s, x))


strat = Strat(k, ["EUR", "AMR"], "ancestry", k.pc)
wts = {s: strat.DEN[:, strat.j(s, D)].sum() for s in ["EUR", "AMR"]}
tw = sum(wts.values())
wts = {s: v / tw for s, v in wts.items()}
ncls = {s: (strat.ncl(s, D), strat.ncl(s, E)) for s in ["EUR", "AMR"]}
log("FAMILY C units", len(strat.units), "weights", wts, "clusters", ncls,
    "sets", {s: (int(((k.ancestry == s) & (k.trait_scope == D)).sum()), int(((k.ancestry == s) & (k.trait_scope == E)).sum()))
             for s in ["EUR", "AMR"]})
ob = strat.obs()


def family_c(g):
    f, den = strat.draws(g)
    res = {}
    defs = {"C2_EUR": {"EUR": 1.0}, "C3_AMR": {"AMR": 1.0}, "C1_EURAMR": wts}
    for name, w in defs.items():
        pt = sum(v * (ob[strat.j(s, D)] - ob[strat.j(s, E)]) for s, v in w.items())
        dr = sum(v * (f[:, strat.j(s, D)] - f[:, strat.j(s, E)]) for s, v in w.items())
        empty = np.zeros(B, bool)
        for s in w:
            for x in (D, E):
                empty |= den[:, strat.j(s, x)] == 0
        ok = np.isfinite(dr)
        fin = dr[ok]
        sd = float(np.std(fin, ddof=1))
        fd = sum(v * ob[strat.j(s, D)] for s, v in w.items())
        fe = sum(v * ob[strat.j(s, E)] for s, v in w.items())
        res[name] = dict(fraction_direct=float(fd), fraction_enzyme=float(fe), point=float(pt),
                         lo95=float(np.percentile(fin, 2.5)), hi95=float(np.percentile(fin, 97.5)),
                         lo98=float(np.percentile(fin, LO98)), hi98=float(np.percentile(fin, HI98)),
                         n_finite=int(ok.sum()), n_empty_arm=int(empty.sum()), sd=sd, mde=MDE_FACTOR * sd)
    return res


famc = {}
for lab, seed in [("same_seed_20260923", 20260923), ("alt_seed_1", 1), ("alt_seed_2", 2)]:
    famc[lab] = family_c(gen(seed, "H58:C:boot"))
    for mname, v in famc[lab].items():
        log("FAMILY C", lab, mname, json.dumps({a: round(b, 6) if isinstance(b, float) else b for a, b in v.items()}))

# non-EUR pooled descriptive row (separate pooled clusters over non-EUR sets)
ne = k[k.ancestry != "EUR"].copy()
ne["pc"] = single_linkage(ne)
cm_ne = class_cluster_means(ne)
fr_ne = {x: cm_ne.loc[cm_ne.trait_scope == x, "num"].sum() / cm_ne.loc[cm_ne.trait_scope == x, "den"].sum() for x in (D, E)}
wide_ne = cm_ne.pivot(index="pc", columns="trait_scope", values=["num", "den"]).fillna(0.0)
wsave = wide
wide = wide_ne
ne_boot = pooled_boot(sorted(wide_ne.index), gen(20260923, "H58:desc:nonEUR_pooled"))
wide = wsave
log(f"NONEUR descriptive {100*fr_ne[D]:.4f} vs {100*fr_ne[E]:.4f} diff {fr_ne[D]-fr_ne[E]:.6f} boot {json.dumps(ne_boot)} "
    f"clusters {cm_ne[cm_ne.trait_scope==D].shape[0]}/{cm_ne[cm_ne.trait_scope==E].shape[0]}")


# ------------------------------------------------------------------ 5. calibration: random splits of one-class clusters
def arm(fr):
    fr = fr.copy()
    fr["pc"] = single_linkage(fr)
    return fr, sorted(fr.pc.unique())


eas, eas_ids = arm(k[(k.ancestry == "EAS") & (k.trait_scope == E)])
afr, afr_ids = arm(k[(k.ancestry == "AFR") & (k.trait_scope == E)])
log("calibration pools: EAS enzyme clusters", len(eas_ids), "AFR enzyme clusters", len(afr_ids))


def split(fr, ids, perm, nd, ne_, stratum):
    dset = {ids[i] for i in perm[:nd]}
    eset = {ids[i] for i in perm[nd:nd + ne_]}
    sub = fr[fr.pc.isin(dset | eset)].copy()
    sub["trait_scope"] = np.where(sub.pc.isin(dset), D, E)
    sub["ps"] = stratum
    return sub.drop(columns="pc")


def covers(fin, lo, hi):
    a, b = np.percentile(fin, lo), np.percentile(fin, hi)
    return bool(a <= 0 <= b)


def calib(seed):
    (d2, e2), (d3, e3) = ncls["EUR"], ncls["AMR"]
    afr_e = min(e3, len(afr_ids) - d3)
    out = {"N1a": [0, 0, 10 ** 9], "N1b": [0, 0, 10 ** 9], "N1c": [0, 0, 10 ** 9]}
    for i in range(NSPLIT):
        fa = split(eas, eas_ids, gen(seed, f"H58:N1a:split:{i}").permutation(len(eas_ids)), d2, e2, "EAS")
        fb = split(eas, eas_ids, gen(seed, f"H58:N1b:split:{i}").permutation(len(eas_ids)), d3, e3, "EAS")
        fc_afr = split(afr, afr_ids, gen(seed, f"H58:N1c:splitAFR:{i}").permutation(len(afr_ids)), d3, afr_e, "AFR")
        fc = pd.concat([fa, fc_afr], ignore_index=True)
        for name, fr, strata, w in [("N1a", fa, ["EAS"], {"EAS": 1.0}), ("N1b", fb, ["EAS"], {"EAS": 1.0}),
                                    ("N1c", fc, ["EAS", "AFR"], {"EAS": wts["EUR"], "AFR": wts["AMR"]})]:
            fr = fr.reset_index(drop=True)
            st = Strat(fr, strata, "ps", single_linkage(fr))
            f, _ = st.draws(gen(seed, f"H58:{name}:boot:{i}"))
            dr = sum(v * (f[:, st.j(s, D)] - f[:, st.j(s, E)]) for s, v in w.items())
            fin = dr[np.isfinite(dr)]
            out[name][0] += covers(fin, 2.5, 97.5)
            out[name][1] += covers(fin, LO98, HI98)
            out[name][2] = min(out[name][2], len(fin))
    return {n: dict(cover95=v[0], cover98=v[1], min_finite=v[2], passes=(v[0] >= 180 and v[1] >= 190)) for n, v in out.items()}


cal = {}
CAL_SEEDS = [(l, s) for l, s in [("same_seed_20260923", 20260923), ("alt_seed_1", 1), ("alt_seed_2", 2)]
             if l in os.environ.get("U5_CAL", "same_seed_20260923,alt_seed_1,alt_seed_2").split(",")]
for lab, seed in (CAL_SEEDS if NSPLIT else []):
    cal[lab] = calib(seed)
    log("CALIBRATION", lab, json.dumps(cal[lab]))


# ------------------------------------------------------------------ 6. coverage check (lead windows with a credible set)
cred = pd.read_csv(CREDSETS, usecols=["study", "locus", "chromosome", "position"], dtype={"locus": str})


def nk(s):
    return s.astype(str).str.strip().str.rstrip("0").str.rstrip(".")


cred["k"] = nk(cred.locus)
cov_rows = []
for (anc, sc, st), _ in k.groupby(["ancestry", "trait_scope", "study"]):
    t = L[(L.study == st) & (L.index_p < GWS)].copy()
    c = cred[cred.study == st]
    keys = set(c.k)
    nw = nkey = npos = 0
    if len(t):
        t["k"] = nk(t.locus)
        for ch, g in t.sort_values(["CHR", "BP"]).groupby("CHR"):
            bp = g.BP.to_numpy()
            wid = np.cumsum(np.concatenate([[True], np.diff(bp) > GAP]))
            cp = c.loc[c.chromosome.astype(str) == str(ch), "position"].to_numpy()
            for w_ in np.unique(wid):
                sel = wid == w_
                nw += 1
                nkey += bool(set(g.k.to_numpy()[sel]) & keys)
                lo, hi = bp[sel].min() - 500_000, bp[sel].max() + 500_000
                npos += bool(((cp >= lo) & (cp <= hi)).any())
    cov_rows.append(dict(ancestry=anc, trait_scope=sc, study=st, windows=nw, with_cs_key=nkey, with_cs_pos=npos))
cov = pd.DataFrame(cov_rows)
cov.to_csv(OUT / "coverage_per_study.tsv", sep="\t", index=False)
covg = cov.groupby(["ancestry", "trait_scope"])[["windows", "with_cs_key", "with_cs_pos"]].sum().reset_index()
covg.to_csv(OUT / "coverage.tsv", sep="\t", index=False)
log("COVERAGE\n" + covg.to_string(index=False))
cov_pass = {}
for s in ["EUR", "AMR"]:
    sh = [float(covg[(covg.ancestry == s) & (covg.trait_scope == x)].eval("with_cs_key / windows").iloc[0]) for x in (D, E)]
    cov_pass[s] = bool(min(sh) > 0 and max(sh) <= 2 * min(sh))
log("COVERAGE pass", cov_pass)


# ------------------------------------------------------------------ 7. calls
def calls(fc, cl):
    res = {}
    for name, cov_ok, ctl in [("C2_EUR", cov_pass["EUR"], "N1a"),
                              ("C3_AMR", cov_pass["AMR"], "N1b"),
                              ("C1_EURAMR", cov_pass["EUR"] and cov_pass["AMR"], "N1c")]:
        v = fc[name]
        reasons = []
        if not cov_ok:
            reasons.append("coverage")
        if v["mde"] > diff:
            reasons.append(f"power (MDE {v['mde']:.4f} > {diff:.4f})")
        if not cl[ctl]["passes"]:
            reasons.append(f"calibration {ctl} {cl[ctl]['cover95']}/{cl[ctl]['cover98']}")
        if v["n_empty_arm"] > 100:
            reasons.append(f"empty arm {v['n_empty_arm']}")
        if reasons:
            res[name] = "no call: " + "; ".join(reasons)
        elif v["lo98"] > 0:
            res[name] = "difference"
        else:
            res[name] = "no difference detected"
    return res


for lab in cal:
    log("CALLS", lab, json.dumps(calls(famc[lab], cal[lab])))

summary = dict(counts=counts, clusters=dict(total=ncl, direct=len(dcl), enzyme=len(ecl), shared=len(dcl & ecl),
                                            direct_with_eur=len(dcl & eur_d_cl)),
               pooled=dict(fraction_direct=frac[D], fraction_enzyme=frac[E], difference=diff), composition=comp,
               pooled_boot=pool_res, family_c_units=len(strat.units), weights=wts, family_c_clusters=ncls,
               family_c=famc, calibration=cal, coverage_pass=cov_pass,
               nonEUR=dict(fd=fr_ne[D], fe=fr_ne[E], boot=ne_boot),
               mde_factor=MDE_FACTOR, B=B, n_split=NSPLIT,
               inputs={str(p): sha(p) for p in [COORDS, ARCH, CREDSETS]})
(OUT / "summary.json").write_text(json.dumps(summary, indent=1, default=float))
(OUT / "log.txt").write_text("\n".join(LOG) + "\n")
