#!/usr/bin/env python
"""U15 independent retabulation of Figure S2A, S2B, S2D (values), S2F (left), S2G, S2H, S2I counts.

Written from the S2 legends (docs/manuscript/05_figure_legends.md) and the producers' input files.
Mode "coloc" reads only COLOC r2 + credible_sets.csv (light); mode "carma" adds the S2F right
half from combined_finemapping.csv (heavy, sbatch only).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
FM = ROOT / "GWAS/finemapping/results"
COLOC = FM / "susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv"
TIER = ROOT / "GWAS/finemapping/config/gwas_trait_tier.tsv"
CRED = FM / "credible_sets.csv"
ANNOT = ROOT / "figures/candidates/fig2-coloc-r2-20260924T053244Z/variant_class/coloc_variant_annotation.csv"
DEG = ROOT / "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv"
OUT = Path(sys.argv[1])
MODE = sys.argv[2] if len(sys.argv) > 2 else "coloc"
OUT.mkdir(parents=True, exist_ok=True)
RES = {}
LOG = []


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.append(s)


tier = pd.read_csv(TIER, sep="\t")
MAIN = sorted(tier.loc[tier.placement == "main", "study_name"])
assert len(MAIN) == 35
DIRECT = set(tier.loc[(tier.placement == "main") & (tier.tier == 1), "study_name"])
log("main strata", len(MAIN), "direct (tier1)", len(DIRECT))

# ------------------------------------------------------------------ COLOC r2, 35 strata
co = pd.read_csv(COLOC, usecols=["gwas_name", "gene", "ensembl", "chr", "PP.H4.abf", "PP.H4.susie", "top_snp",
                                 "ancestry", "method"], dtype={"chr": str})
n_all = len(co)
co = co[co.gwas_name.isin(MAIN)].copy()
log("COLOC rows all", n_all, "in 35 strata", len(co), "strata present", co.gwas_name.nunique())
named = co.gene.notna() & (co.gene.astype(str).str.strip() != "")
log("rows without a gene name", int((~named).sum()))
dup = co[named].duplicated(["gwas_name", "gene"]).sum()
log("duplicated (stratum, gene name) rows", int(dup), "duplicated (stratum, ensembl)", int(co.duplicated(["gwas_name", "ensembl"]).sum()))
anc_by_study = co.groupby("gwas_name").ancestry.agg(lambda s: sorted(s.dropna().unique()))
log("ancestry labels per stratum:", {k: v for k, v in anc_by_study.items()})
co["anc"] = co.gwas_name.map(lambda s: anc_by_study[s][0] if len(anc_by_study[s]) == 1 else "MIXED")
co.loc[co.anc == "CSA", "anc"] = "SAS"
su = co["PP.H4.susie"]
ab = co["PP.H4.abf"]
cn = co[named]

# S2A right: genes per stratum with either PP.H4 > 0.5
pos = cn[(cn["PP.H4.susie"] > 0.5) | (cn["PP.H4.abf"] > 0.5)]
per = pos.groupby("gwas_name").gene.nunique()
svein = [s for s in MAIN if s.startswith("2023_36280732_NAFLD_")]
rows = per.drop(svein, errors="ignore").to_dict()
rows["Sveinbjornsson (3 cohorts)"] = pos[pos.gwas_name.isin(svein)].gene.nunique()
RES["S2A_genes_per_row"] = rows
RES["S2A_sum_rows"] = int(sum(rows.values()))
RES["S2A_union"] = int(pos.gene.nunique())
RES["S2A_strata_with_no_positive"] = [s for s in MAIN if s not in per.index]
log("S2A rows", len(rows), "sum", RES["S2A_sum_rows"], "union", RES["S2A_union"])

# S2F left
both = co[su.notna() & ab.notna()]
RES["S2F_pairs_both"] = int(len(both))
RES["S2F_r"] = float(np.corrcoef(both["PP.H4.abf"], both["PP.H4.susie"])[0, 1])
RES["S2F_abf_pairs_total"] = int(ab.notna().sum())
RES["S2F_susie_pairs_total"] = int(su.notna().sum())
RES["S2F_pairs_susie_gt05"] = int((su > 0.5).sum())
RES["S2F_pairs_abf_gt05"] = int((ab > 0.5).sum())
RES["S2F_genes_susie_gt05"] = int(cn.loc[cn["PP.H4.susie"] > 0.5, "gene"].nunique())
RES["S2F_genes_abf_gt05"] = int(cn.loc[cn["PP.H4.abf"] > 0.5, "gene"].nunique())
RES["S2F_genes_susie_gt05_ensembl"] = int(co.loc[su > 0.5, "ensembl"].nunique())
log("S2F", {k: v for k, v in RES.items() if k.startswith("S2F")})

# S2G: genes by threshold x method, EUR if any EUR stratum passes, else non-EUR only
g_rows = []
for t in (0.5, 0.8, 0.9):
    for meth, col in (("ABF", "PP.H4.abf"), ("SuSiE", "PP.H4.susie")):
        p = cn[cn[col] > t]
        genes = set(p.gene)
        eur = set(p.loc[p.anc == "EUR", "gene"])
        g_rows.append(dict(threshold=t, method=meth, total=len(genes), EUR=len(eur), nonEUR_only=len(genes - eur)))
RES["S2G"] = g_rows
log("S2G\n" + pd.DataFrame(g_rows).to_string(index=False))

# S2H / S2I
cs_ = cn[cn["PP.H4.susie"].notna()].copy()
ba = cs_.groupby(["gene", "anc"])["PP.H4.susie"].max().unstack()
best_idx = cs_.groupby("gene")["PP.H4.susie"].idxmax()
best = cs_.loc[best_idx, ["gene", "ensembl", "chr", "top_snp", "gwas_name", "PP.H4.susie"]].set_index("gene")
best["n_anc09"] = (ba > 0.9).sum(axis=1)
best["n_anc05"] = (ba > 0.5).sum(axis=1)
hi = best[best["PP.H4.susie"] > 0.9]
RES["S2H_genes_gt09"] = int(len(hi))
RES["S2H_chr_counts"] = hi.chr.value_counts().sort_index().to_dict()
RES["S2H_labeled_n_anc09_ge2"] = int((hi.n_anc09 >= 2).sum())
s05 = best[best["PP.H4.susie"] > 0.5]
RES["S2I_genes_gt05"] = int(len(s05))
# ties for best row: report how many genes have >1 row at the max
nties = cs_.merge(best[["PP.H4.susie"]].rename(columns={"PP.H4.susie": "mx"}), left_on="gene", right_index=True)
nties = nties[nties["PP.H4.susie"] == nties.mx].groupby("gene").size()
RES["S2HI_genes_with_tied_best_rows"] = int((nties > 1).sum())

# lead-variant class from the variant annotation used by the producer
an = pd.read_csv(ANNOT)
log("annotation columns", list(an.columns), "rows", len(an))
RES["annot_columns"] = list(an.columns)
an.to_csv(OUT / "annot_head.tsv", sep="\t", index=False) if len(an) < 5000 else an.head(50).to_csv(OUT / "annot_head.tsv", sep="\t", index=False)
best.reset_index().to_csv(OUT / "S2HI_best_per_gene.tsv", sep="\t", index=False)
ba.reset_index().to_csv(OUT / "S2HI_best_per_gene_ancestry.tsv", sep="\t", index=False)

# bulk calls joined on Ensembl ID (version stripped), gene universe 24,115
deg = pd.read_csv(DEG)
deg["ens"] = deg.gene.str.replace(r"\.\d+$", "", regex=True)
assert deg.ens.is_unique
RES["bulk_universe"] = int(len(deg))
hh = hi.reset_index().merge(deg[["ens", "logFC", "padj"]], left_on="ensembl", right_on="ens", how="left")
call = np.where(hh.padj.isna(), "not tested",
                np.where((hh.padj < 0.05) & (hh.logFC.abs() > 0.5), np.where(hh.logFC > 0, "Up", "Down"), "n.s."))
RES["S2H_bulk_calls"] = pd.Series(call).value_counts().to_dict()
log("S2H", RES["S2H_genes_gt09"], RES["S2H_labeled_n_anc09_ge2"], RES["S2H_bulk_calls"], "chr", RES["S2H_chr_counts"])

# S2D values
ef = co[(co.gene == "EFHD1") & co.gwas_name.isin(["UKBB_AST", "BBJ_AST"])][["gwas_name", "PP.H4.susie", "PP.H4.abf", "top_snp"]]
RES["S2D_EFHD1"] = ef.to_dict("records")
d_ef = deg[deg.symbol == "EFHD1"][["gene", "logFC", "padj"]]
RES["S2D_EFHD1_bulk"] = d_ef.to_dict("records")

# ------------------------------------------------------------------ S2A left, S2B: credible_sets.csv
cr = pd.read_csv(CRED, usecols=["chromosome", "position", "locus", "study", "susie_pip", "carma_pip"], low_memory=False)
cr = cr[cr.study.isin(MAIN)]
loci = cr.groupby("study").locus.nunique()
RES["S2A_loci_per_stratum"] = loci.to_dict()
RES["S2A_svein_loci_union"] = int(cr[cr.study.isin(svein)].locus.nunique())
RES["S2A_fine_mapped_strata"] = int(loci.size)
win = cr.groupby(["study", "locus"]).agg(chr=("chromosome", "first"), lo=("position", "min"), hi=("position", "max"),
                                         nchr=("chromosome", "nunique")).reset_index()
assert (win.nchr == 1).all()
RES["S2B_windows"] = int(len(win))


def merge_windows(w, pad):
    out = []
    for ch, g in w.assign(a=w.lo - pad, b=w.hi + pad).sort_values(["chr", "a"]).groupby("chr"):
        cur = None
        for r in g.itertuples(index=False):
            if cur is None or r.a > cur["b"]:
                if cur is not None:
                    out.append(cur)
                cur = dict(chr=ch, a=r.a, b=r.b, strata={r.study}, keys={(r.study, r.locus)})
            else:
                cur["b"] = max(cur["b"], r.b)
                cur["strata"].add(r.study)
                cur["keys"].add((r.study, r.locus))
        out.append(cur)
    return out


m250 = merge_windows(win, 250_000)
m0 = merge_windows(win, 0)
RES["S2B_merged_250kb"] = len(m250)
RES["S2B_one_stratum"] = sum(len(x["strata"]) == 1 for x in m250)
RES["S2B_two_plus"] = sum(len(x["strata"]) >= 2 for x in m250)
RES["S2B_merged_nopad"] = len(m0)
key2m = {k: i for i, x in enumerate(m250) for k in x["keys"]}
cr["m"] = [key2m[(s, l)] for s, l in zip(cr.study, cr.locus)]
bs = cr.groupby("m").agg(susie=("susie_pip", "max"), carma=("carma_pip", "max"))


def bands(x):
    x = x.dropna()
    return dict(ge09=int((x >= 0.9).sum()), b05_09=int(((x >= 0.5) & (x < 0.9)).sum()), lt05=int((x < 0.5).sum()), n=int(len(x)))


RES["S2B_susie"] = bands(bs.susie)
RES["S2B_carma"] = bands(bs.carma)
sx = pd.read_csv(FM / "susiex_mvp/susiex_cs_mvp.csv")
sx = sx[sx.trait_pair.isin(["ALT", "AST", "NAFLD"])]
RES["S2B_susiex"] = bands(sx.groupby(["locus_id", "trait_pair"]).OVRL_PIP.max())
me = pd.read_csv(FM / "mesusie_mvp/mesusie_locus_summary_mvp.csv")
me_all = me[me.trait_pair.isin(["ALT", "AST", "NAFLD"])]
RES["S2B_mesusie_converged"] = bands(me_all[me_all.converged == True].max_pip)  # noqa: E712
RES["S2B_mesusie_all"] = bands(me_all.max_pip)
log("S2B", {k: v for k, v in RES.items() if k.startswith("S2B")})
# EFHD1 UKBB_AST credible set near chr2:233.5 Mb (untargeted run)
RES["S2D_UKBB_AST_rows_chr2_233Mb"] = int(((cr.study == "UKBB_AST") & (cr.chromosome.astype(str) == "2")
                                          & cr.position.between(233_300_000, 233_800_000)).sum())
RES["S2D_BBJ_AST_rows_chr2_233Mb"] = int(((cr.study == "BBJ_AST") & (cr.chromosome.astype(str) == "2")
                                         & cr.position.between(233_300_000, 233_800_000)).sum())

# ------------------------------------------------------------------ S2F right (heavy)
if MODE == "carma":
    n_all = n_35 = 0
    parts = []
    for ch in pd.read_csv(FM / "combined_finemapping.csv", usecols=["study", "susie_pip", "carma_pip", "both_in_cs",
                                                                     "either_in_cs", "carma_outlier"],
                          dtype={"both_in_cs": str, "either_in_cs": str, "carma_outlier": str}, chunksize=2_000_000):
        n_all += len(ch)
        ch = ch[ch.study.isin(MAIN)]
        n_35 += len(ch)
        parts.append(ch[ch.susie_pip.notna() & ch.carma_pip.notna()])
    RES["S2F_cf_rows_all"] = int(n_all)
    RES["S2F_cf_rows_35"] = int(n_35)
    b2 = pd.concat(parts, ignore_index=True)
    RES["S2F_variants_both"] = int(len(b2))
    RES["S2F_carma_strata"] = sorted(b2.study.unique())
    RES["S2F_r_susie_carma"] = float(np.corrcoef(b2.susie_pip, b2.carma_pip)[0, 1])

    def tf(s):
        return s.astype(str).str.upper().isin(["TRUE", "1", "1.0"])
    RES["S2F_bothCS_over_eitherCS"] = float(tf(b2.both_in_cs).sum() / tf(b2.either_in_cs).sum())
    RES["S2F_carma_outlier_frac"] = float(tf(b2.carma_outlier).mean())
    log("S2F right", {k: v for k, v in RES.items() if k.startswith("S2F_") and ("cf" in k or "carma" in k or "both" in k)})

(OUT / f"results_{MODE}.json").write_text(json.dumps(RES, indent=1, default=lambda o: int(o) if isinstance(o, np.integer) else str(o)))
(OUT / f"log_{MODE}.txt").write_text("\n".join(LOG) + "\n")
