"""Independent check of claim R02-03 (Fig 2A cascade loci vs COLOC genes).

Reads only saved tables. Reimplements the Fig 2A loci logic from
scripts/figures/fig1_gwas_creative_options.py (_unique_loci_by_at, _load_manifest)
without importing it, then asks where the 436 multi-signal genes' positive
COLOC leads sit relative to fine-mapped credible-set windows.
"""
import csv
from collections import defaultdict
import pandas as pd

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
TIER = f"{ROOT}/GWAS/finemapping/config/gwas_trait_tier.tsv"
REG = f"{ROOT}/GWAS/finemapping/config/gwas_registry.tsv"
CS = f"{ROOT}/GWAS/finemapping/results/credible_sets.csv"
COLOC = f"{ROOT}/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv"
ARCH = (f"{ROOT}/GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/"
        "run_19783321/credible_set_pip_architecture.tsv")
PAD = 250_000

TOK = [("ChronLiver", "CLD"), ("Cirrhosis", "Cirr"), ("Albumin", "Alb"), ("Platelet", "Plt"),
       ("PDFF", "PDFF"), ("NAFLD", "NAFLD"), ("NASH", "NASH"), ("HCC", "HCC"),
       ("ALT", "ALT"), ("AST", "AST"), ("GGT", "GGT")]


def trait_of(s):
    lo = s.lower()
    for t, lab in TOK:
        if t.lower() in lo:
            return lab
    return None


MAIN = {r["study_name"] for r in csv.DictReader(open(TIER), delimiter="\t") if r["placement"] == "main"}
anc = {r["study_name"]: r["ancestry"] for r in csv.DictReader(open(REG), delimiter="\t")
       if r["study_name"] in MAIN}
print("MAIN strata:", len(MAIN))

# ---- credible-set windows per (study, locus), MAIN only (same as producer) ----
cs = pd.read_csv(CS, dtype={"chromosome": str, "locus": str, "study": str})
cs = cs[cs.study.isin(MAIN)].copy()
cs["pos"] = cs.position.astype(float).astype(int)
cs["pip"] = pd.to_numeric(cs.recommended_pip, errors="coerce").fillna(
    pd.to_numeric(cs.max_pip, errors="coerce")).fillna(0.0)
win = cs.groupby(["study", "locus"]).agg(chrom=("chromosome", "first"), mn=("pos", "min"),
                                         mx=("pos", "max"), pip=("pip", "max")).reset_index()
win["anc"] = win.study.map(anc)
win["trait"] = win.study.map(trait_of)
print("MAIN studies with credible sets:", win.study.nunique())
print("MAIN studies WITHOUT credible sets:", sorted(MAIN - set(win.study)))
print("raw (study,locus) windows per UKBB enzyme study:",
      win[win.study.str.startswith("UKBB_")].groupby("study").size().to_dict())


def merge(ivs, pad):
    lst = sorted((mn - pad, mx + pad, pip, a) for mn, mx, pip, a in ivs)
    out = []
    for ps, pe, pip, a in lst:
        if out and ps <= out[-1][0]:
            if pe > out[-1][0]:
                out[-1][0] = pe
            if pip > out[-1][1]:
                out[-1][1], out[-1][2] = pip, a
        else:
            out.append([pe, pip, a])
    return [(m[1], m[2]) for m in out]


at = defaultdict(int)
for (t, c), g in win.groupby(["trait", "chrom"]):
    for _p, a in merge(list(zip(g.mn, g.mx, g.pip, g.anc)), PAD):
        at[(a, t)] += 1
glob = sum(len(merge(list(zip(g.mn, g.mx, g.pip, g.anc)), PAD)) for _c, g in win.groupby("chrom"))
print("within-trait loci total:", sum(at.values()), " global unique:", glob)
for t in ["NAFLD", "NASH", "PDFF", "ALT", "AST", "GGT"]:
    print("  ", t, sum(v for (a, tt), v in at.items() if tt == t),
          {a: v for (a, tt), v in at.items() if tt == t})

# ---- multi-signal genes (max PP.H4.susie > 0.5 over MAIN rows, blank gene excluded) ----
co = pd.read_csv(COLOC, usecols=["gwas_name", "gene", "chr", "PP.H4.susie", "PP.H4.abf",
                                 "top_snp", "method"], dtype={"chr": str, "gene": str})
co = co[co.gwas_name.isin(MAIN)].copy()
co["gene"] = co.gene.fillna("")
co = co[co.gene.str.strip() != ""]
co["su"] = pd.to_numeric(co["PP.H4.susie"], errors="coerce")
co["abf"] = pd.to_numeric(co["PP.H4.abf"], errors="coerce")
gmax_su = co.groupby("gene").su.max()
gmax_abf = co.groupby("gene").abf.max()
ms = set(gmax_su[gmax_su > 0.5].index)
union = ms | set(gmax_abf[gmax_abf > 0.5].index)
print("multi-signal genes:", len(ms), " union:", len(union), " single-only:", len(union - ms))

pos = co[(co.su > 0.5) & co.gene.isin(ms)].copy()
print("positive (gene, study) pairs:", len(pos))
sp = pos.top_snp.astype(str).str.split(":", expand=True)
pos["lchr"] = sp[0].str.replace("chr", "", regex=False)
pos["lpos"] = pd.to_numeric(sp[1], errors="coerce")
print("positive pairs with unparseable top_snp:", int(pos.lpos.isna().sum()))

# build window index
wi = defaultdict(list)          # (study, chrom) -> [(mn,mx)]
wi_any = defaultdict(list)      # chrom -> [(mn,mx)]
for r in win.itertuples():
    wi[(r.study, r.chrom)].append((r.mn, r.mx))
    wi_any[r.chrom].append((r.mn, r.mx))


def inside(lst, p, pad):
    return any(mn - pad <= p <= mx + pad for mn, mx in lst)


for pad in (0, PAD):
    same, anyw = set(), set()
    for r in pos.itertuples():
        if pd.isna(r.lpos):
            continue
        if inside(wi.get((r.gwas_name, r.lchr), []), r.lpos, pad):
            same.add(r.gene)
        if inside(wi_any.get(r.lchr, []), r.lpos, pad):
            anyw.add(r.gene)
    print(f"pad={pad}: genes with a positive lead in a same-study window: {len(same)}; "
          f"in any MAIN window: {len(anyw)}")

MVP_EUR = {"MVP_ALT_EUR", "MVP_AST_EUR", "MVP_NAFLD_EUR"}
UKBB_ENZ = {"UKBB_ALT", "UKBB_AST", "UKBB_GGT"}
pos_sets = pos.groupby("gene").gwas_name.agg(set)
only_mvp = [g for g, s in pos_sets.items() if s <= MVP_EUR]
only_mvp_ukbb = [g for g, s in pos_sets.items() if s <= (MVP_EUR | UKBB_ENZ)]
any_mvp = [g for g, s in pos_sets.items() if s & MVP_EUR]
no_cs_study = [g for g, s in pos_sets.items() if not (s & set(win.study))]
print("genes positive ONLY in MVP EUR strata:", len(only_mvp))
print("genes positive ONLY in MVP EUR + UKBB enzyme strata:", len(only_mvp_ukbb))
print("genes positive in >=1 MVP EUR stratum:", len(any_mvp))
print("genes with no positive study that has any credible set:", len(no_cs_study))

# strongest-SuSiE study per gene (producer's assignment), tie -> smaller gwas_name
best = (pos.sort_values(["gene", "su", "gwas_name"], ascending=[True, False, True])
        .drop_duplicates("gene"))
print("strongest-study counts, top:", best.gwas_name.value_counts().head(12).to_dict())
eur_alt = best[best.gwas_name.isin(["MVP_ALT_EUR", "UKBB_ALT"])]
print("EUR ALT multi-signal genes by strongest study:", eur_alt.gwas_name.value_counts().to_dict())
eur_enz = best[best.gwas_name.map(lambda s: anc.get(s) == "EUR" and trait_of(s) in ("ALT", "AST", "GGT"))]
print("EUR enzyme multi-signal genes (strongest study):", len(eur_enz),
      eur_enz.gwas_name.value_counts().to_dict())

# ---- Fig 2D arm composition ----
ar = pd.read_csv(ARCH, sep="\t")
ar_pass = ar[ar.pip_sum_gate_passed.astype(str).str.lower() == "true"]
print("2D instances total / passing:", len(ar), len(ar_pass))
print("2D studies per scope:", ar_pass.groupby("trait_scope").study.nunique().to_dict())
main_direct = {s for s in MAIN if trait_of(s) in ("NAFLD", "NASH", "PDFF")}
main_enz = MAIN - main_direct
print("MAIN direct strata absent from 2D direct arm:",
      sorted(main_direct - set(ar_pass[ar_pass.trait_scope == "tier1_direct_masld_pdff"].study)))
print("MAIN enzyme strata absent from 2D enzyme arm:",
      sorted(main_enz - set(ar_pass[ar_pass.trait_scope == "tier2_liver_enzyme"].study)))
eur_enz_arm = ar_pass[(ar_pass.trait_scope == "tier2_liver_enzyme") & (ar_pass.ancestry == "EUR")]
print("2D EUR enzyme-arm passing instances:", len(eur_enz_arm), eur_enz_arm.study.value_counts().to_dict())
eur_enz_all = ar[(ar.trait_scope == "tier2_liver_enzyme") & (ar.ancestry == "EUR")]
print("2D EUR enzyme-arm all instances:", len(eur_enz_all))
