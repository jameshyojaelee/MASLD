"""DIG1 job 41b (exploratory; follows job 40). Genome-wide: for genes whose liver eQTL colocalizes with a
liver-injury trait, does established disease move liver expression in the direction that the genetics
links to more injury (concordant) or the opposite direction (discordant)? The Resource compares no
directions; this is a new exploratory question.
Genetic direction per pair: sign(z_eQTL x z_GWAS) at the colocalized SNP, alleles aligned (eQTL EA/NEA
vs GWAS allele1 = effect allele; palindromic SNPs dropped). Traits: tier 1 (NAFLD, NASH, PDFF) and
tier 2 (ALT, AST, GGT) of gwas_trait_tier.tsv; higher = worse for all. +1 = the allele that raises
liver expression raises the trait. A sign needs |z_eQTL| >= 3 and |z_GWAS| >= 2 at the SNP.
Pairs: SuSiE signal pairs PP.H4 > 0.8 (SNP = eQTL signal hit2); ABF pairs PP.H4.abf > 0.8 (SNP = top_snp)
as a sensitivity. Gene direction = the pair with the highest PP4 (primary) or the majority sign (mixed
genes dropped).
Disease direction: job 41a (pooled disease-vs-control logFC of the adopted 774-unit fit; NAS and
fibrosis partial Spearman within disease units).
Independence: genes are clustered by colocalized SNP position (greedy by PP4, 1 Mb on the same
chromosome); one gene per cluster.
Primary family (Holm over 2), Fisher exact test of genetic sign vs disease sign, SuSiE pairs, all
ancestries, top-PP4 rule, one gene per locus:
  K1  genes with disease-vs-control padj < 0.05, sign of logFC.
  K2  genes with NAS BH < 0.05, sign of rho.
Descriptive: the same over variants (ABF, union, EUR-only, majority rule, no dedup, enzyme-only,
MASLD-trait-only, called genes |logFC| > 0.5, fibrosis, all genes by sign of t). Positive control:
z of the PNPLA3 rs738409 G allele (hg19 22:44324727) and TM6SF2 rs58542926 T allele (19:19379549) in
every trait [hg19 positions recalled]. Seed not needed (no resampling).
Usage: python -I 41b_genetic_vs_disease_direction.py <new_out_dir>
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

# eQTL lookups
need = P[["chr", "ensembl", "pos"]].drop_duplicates()
eqr = []
for c, nd in need.groupby("chr"):
    keys = set(zip(nd.ensembl, nd.pos))
    for ch in pd.read_csv(EQ / f"chr{c}_marginal_summary_results.tsv", sep="\t", chunksize=3_000_000,
                          usecols=["CHR", "POS", "NEA", "EA", "Beta", "SE", "PVAL", "ENSG"]):
        ch = ch[ch.POS.isin(nd.pos)]
        ch = ch[[(e, p) in keys for e, p in zip(ch.ENSG, ch.POS)]]
        eqr.append(ch)
    print("eqtl chr", c, flush=True)
eqr = pd.concat(eqr).rename(columns={"ENSG": "ensembl", "POS": "pos", "CHR": "chr"})
eqr["z_eqtl"] = eqr.Beta / eqr.SE
# GWAS lookups
gwr, ctrl = [], []
for s, ps in P.groupby("gwas"):
    pos = set(zip(ps.chr, ps.pos)) | {(c, p) for c, p, _ in CTRL.values()}
    chrs = {c for c, _ in pos}; posset = {p for _, p in pos}
    for ch in pd.read_csv(FM / reg.loc[s, "sumstats_path"], sep="\t", chunksize=4_000_000, dtype={"chromosome": str}):
        ch = ch[ch.position.isin(posset)]
        ch = ch[[(int(c) if str(c).isdigit() else -1, p) in pos for c, p in zip(ch.chromosome, ch.position)]]
        ch = ch.assign(gwas=s)
        gwr.append(ch)
    print("gwas", s, flush=True)
gwr = pd.concat(gwr)
gwr["chr"] = pd.to_numeric(gwr.chromosome, errors="coerce").astype("Int64")
gwr["z_gwas_a1"] = gwr.beta / gwr.se
for name, (c, p, risk) in CTRL.items():
    g = gwr[(gwr.chr == c) & (gwr.position == p)]
    for _, r in g.iterrows():
        z = r.z_gwas_a1 if r.allele1.upper() == risk else (-r.z_gwas_a1 if r.allele2.upper() == risk else np.nan)
        ctrl.append(dict(control=name, gwas=r.gwas, trait=tier.loc[r.gwas, "trait"], z_risk_allele=z, allele1=r.allele1, allele2=r.allele2))
ctrl = pd.DataFrame(ctrl)
ctrl.to_csv(OUT / "positive_control_alleles.tsv", sep="\t", index=False)
bad = set(ctrl[(ctrl.control == "PNPLA3_rs738409_G") & (ctrl.z_risk_allele < 0)].gwas)

M = P.merge(eqr[["chr", "ensembl", "pos", "EA", "NEA", "z_eqtl"]], on=["chr", "ensembl", "pos"], how="left")
M = M.merge(gwr[["gwas", "chr", "position", "allele1", "allele2", "z_gwas_a1"]].rename(columns={"position": "pos"}).astype({"chr": "float"}).assign(chr=lambda d: d.chr.astype("Int64")),
            on=["gwas", "chr", "pos"], how="left")
EA, NEA = M.EA.str.upper(), M.NEA.str.upper(); A1, A2 = M.allele1.str.upper(), M.allele2.str.upper()
same = (EA == A1) & (NEA == A2); swap = (EA == A2) & (NEA == A1)
M["z_gwas"] = np.where(same, M.z_gwas_a1, np.where(swap, -M.z_gwas_a1, np.nan))
M["palindromic"] = [(a, b) in PAL for a, b in zip(EA, NEA)]
M.loc[M.palindromic, "z_gwas"] = np.nan
M["usable"] = M.z_eqtl.abs().ge(3) & M.z_gwas.abs().ge(2)
M["gen_sign"] = np.where(M.usable, np.sign(M.z_eqtl * M.z_gwas), np.nan)
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
