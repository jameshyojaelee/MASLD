"""DIG1 job 41c (post hoc deviation from 41b). 41b excluded every study in which the PNPLA3 rs738409 G allele
had z < 0. In the three GGT studies PNPLA3 has no effect (z -0.14 to -0.87), so they were dropped although
their TM6SF2 control is in the expected direction (UKBB_GGT z +4.1). Here a study fails the control only
if a control risk allele has z < -3. Everything else is as 41b, using 41b's saved pair lookups.
Usage: python -I 41c_direction_posthoc_controls.py <new_out_dir>
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
RUN = ROOT / "Analysis/Hypotheses/DIG1-20261024T113951Z"
FM = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping")
EQ = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/broadaway_eqtl")
CO = ROOT / "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z"
assert len(sys.argv) == 2
OUT = Path(sys.argv[1])
OUT.mkdir(parents=True, exist_ok=False)
PAL = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}
CTRL = {"PNPLA3_rs738409_G": (22, 44324727, "G"), "TM6SF2_rs58542926_T": (19, 19379549, "T")}

tier = pd.read_csv(FM / "config/gwas_trait_tier.tsv", sep="\t").set_index("study_name")
reg = pd.read_csv(FM / "config/gwas_registry.tsv", sep="\t").set_index("study_name")
keep_studies = tier.index[tier.tier <= 2]
sp = pd.read_csv(CO / "susie_coloc_signal_pairs.tsv.gz", sep="\t")
sp = sp[(sp["PP.H4"] > 0.8) & sp.gwas.isin(keep_studies)]
sp = pd.DataFrame(dict(source="susie", gene=sp.gene, ensembl=sp.ensembl, chr=sp.chr.astype(int), gwas=sp.gwas, PP4=sp["PP.H4"],
                       pos=sp.hit2.astype(str).str.split(":").str[1].astype(int)))
ab = pd.read_csv(CO / "susie_coloc_all_gwas.csv", usecols=["gwas_name", "gene", "ensembl", "chr", "PP.H4.abf", "top_snp"])
ab = ab[(pd.to_numeric(ab["PP.H4.abf"], errors="coerce") > 0.8) & ab.gwas_name.isin(keep_studies) & ab.top_snp.notna()]
ab = pd.DataFrame(dict(source="abf", gene=ab.gene, ensembl=ab.ensembl, chr=ab.chr.astype(int), gwas=ab.gwas_name, PP4=ab["PP.H4.abf"].astype(float),
                       pos=ab.top_snp.astype(str).str.split(":").str[1].astype(int)))
P = pd.concat([sp, ab], ignore_index=True)
P["ensembl"] = P.ensembl.astype(str).str.replace(r"\..*$", "", regex=True)
P["ancestry"] = P.gwas.map(reg.ancestry); P["trait"] = P.gwas.map(tier.trait); P["tier"] = P.gwas.map(tier.tier)
print("pairs", P.groupby("source").size().to_dict(), flush=True)


M0 = pd.read_csv(RUN / "out_genetic_vs_disease_direction/pair_directions.tsv.gz", sep="\t")
ctrl = pd.read_csv(RUN / "out_genetic_vs_disease_direction/positive_control_alleles.tsv", sep="\t")
bad = set(ctrl[ctrl.z_risk_allele < -3].gwas)
ctrl.to_csv(OUT / "positive_control_alleles.tsv", sep="\t", index=False)
M = M0.copy()
M["control_failed_study"] = M.gwas.isin(bad)
M.to_csv(OUT / "pair_directions.tsv.gz", sep="\t", index=False)
D = pd.read_csv(RUN / "out_genome_disease_state/gene_disease_state.tsv.gz", sep="\t")


def gene_table(pairs, rule):
    q = pairs[pairs.usable & ~pairs.control_failed_study].sort_values("PP4", ascending=False)
    if rule == "top":
        g = q.drop_duplicates("ensembl")[["ensembl", "gene", "chr", "pos", "PP4", "gwas", "gen_sign"]]
    else:
        agg = q.groupby("ensembl").agg(gene=("gene", "first"), chr=("chr", "first"), pos=("pos", "first"), PP4=("PP4", "max"),
                                       gwas=("gwas", "first"), npos=("gen_sign", lambda s: (s > 0).sum()), nneg=("gen_sign", lambda s: (s < 0).sum()))
        agg = agg[(agg.npos == 0) | (agg.nneg == 0)].reset_index()
        agg["gen_sign"] = np.where(agg.npos > 0, 1.0, -1.0)
        g = agg[["ensembl", "gene", "chr", "pos", "PP4", "gwas", "gen_sign"]]
    return g


def dedup(g):
    g = g.sort_values("PP4", ascending=False); leads = []
    keep = []
    for _, r in g.iterrows():
        if any(r.chr == c and abs(r.pos - p) <= 1_000_000 for c, p in leads):
            keep.append(False)
        else:
            leads.append((r.chr, r.pos)); keep.append(True)
    return g[keep]


def test(g, dcol, filt):
    x = g.merge(D, on="ensembl", how="inner")
    x = x[filt(x)]
    ds = np.sign(x[dcol])
    a = int(((x.gen_sign > 0) & (ds > 0)).sum()); b = int(((x.gen_sign > 0) & (ds < 0)).sum())
    c = int(((x.gen_sign < 0) & (ds > 0)).sum()); d = int(((x.gen_sign < 0) & (ds < 0)).sum())
    n = a + b + c + d
    if n == 0:
        return dict(n=0)
    orr, p = stats.fisher_exact([[a, b], [c, d]])
    pg, pd_ = (a + b) / n, (a + c) / n
    return dict(n=n, gen_up_dis_up=a, gen_up_dis_down=b, gen_down_dis_up=c, gen_down_dis_down=d, frac_concordant=(a + d) / n,
                expected_concordant=pg * pd_ + (1 - pg) * (1 - pd_), odds_ratio=orr, p_fisher=p,
                genes_concordant=",".join(sorted(x.gene[(x.gen_sign * ds) > 0].astype(str))), genes_discordant=",".join(sorted(x.gene[(x.gen_sign * ds) < 0].astype(str))))


DSETS = {"padj<0.05 (logFC)": ("logFC", lambda x: x.padj < 0.05),
         "called |logFC|>0.5 (logFC)": ("logFC", lambda x: (x.padj < 0.05) & (x.logFC.abs() > 0.5)),
         "NAS BH<0.05 (rho)": ("rho_nas", lambda x: x.bh_nas < 0.05),
         "fibrosis BH<0.05 (rho)": ("rho_fib", lambda x: x.bh_fib < 0.05),
         "all genes (sign t)": ("t", lambda x: x.t.notna())}
PSETS = {"susie": M[M.source == "susie"], "abf": M[M.source == "abf"], "union": M}
rows = []
for pn, pm in PSETS.items():
    for anc in ["all", "EUR"]:
        for tr in ["tier1+2", "enzyme_only", "MASLD_only"]:
            q = pm if anc == "all" else pm[pm.ancestry == "EUR"]
            q = q if tr == "tier1+2" else q[q.tier == (2 if tr == "enzyme_only" else 1)]
            for rule in ["top", "majority"]:
                g = gene_table(q, rule)
                for dd in ["locus_dedup", "all_genes"]:
                    gg = dedup(g) if dd == "locus_dedup" else g
                    for dn, (dcol, f) in DSETS.items():
                        rows.append(dict(pairs=pn, ancestry=anc, traits=tr, rule=rule, dedup=dd, disease=dn, n_genes_with_sign=len(gg), **test(gg, dcol, f)))
R = pd.DataFrame(rows)
R.to_csv(OUT / "concordance_grid.tsv", sep="\t", index=False)
base = R[(R.pairs == "susie") & (R.ancestry == "all") & (R.traits == "tier1+2") & (R.rule == "top") & (R.dedup == "locus_dedup")].set_index("disease")
prim = pd.DataFrame([dict(id="K1", disease="padj<0.05 (logFC)", **base.loc["padj<0.05 (logFC)"][["n", "frac_concordant", "expected_concordant", "odds_ratio", "p_fisher"]].to_dict()),
                     dict(id="K2", disease="NAS BH<0.05 (rho)", **base.loc["NAS BH<0.05 (rho)"][["n", "frac_concordant", "expected_concordant", "odds_ratio", "p_fisher"]].to_dict())])
o = np.argsort(prim.p_fisher.to_numpy()); ph = np.maximum.accumulate(prim.p_fisher.to_numpy()[o] * (len(o) - np.arange(len(o))))
prim.loc[prim.index[o], "p_holm"] = np.minimum(1, ph)
prim.to_csv(OUT / "primary_holm.tsv", sep="\t", index=False)
gt = dedup(gene_table(M[M.source == "susie"], "top")).merge(D, on="ensembl", how="left")
gt.to_csv(OUT / "gene_directions_susie_top_dedup.tsv", sep="\t", index=False)
with open(OUT / "pip_freeze.txt", "w") as f:
    import importlib.metadata as im
    f.write("\n".join(sorted(f"{d.metadata['Name']}=={d.version}" for d in im.distributions())))
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 400); pd.set_option("display.max_colwidth", 60)
print(ctrl.to_string())
print("pairs usable", M.groupby("source").usable.mean().round(3).to_dict(), "control-failed studies", sorted(bad))
print(prim.round(4).to_string())
cols = ["pairs", "ancestry", "traits", "rule", "dedup", "disease", "n", "frac_concordant", "expected_concordant", "odds_ratio", "p_fisher"]
print(R[cols].round(3).to_string())
