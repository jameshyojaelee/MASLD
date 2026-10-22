"""Count n for lens-L09 candidates from files on disk. Read-only. Usage: python ncount.py C1|C2|C3|C4|C5|C6"""
import sys
import numpy as np
import pandas as pd
from scipy import stats

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
POOLED = f"{ROOT}/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv"
COHORT = f"{ROOT}/figures/candidates/bulk-context-774-20261014T140904Z_22236311/cohort_disease_all_gene_results.tsv"
ROSTER = f"{ROOT}/figures/candidates/resource-corrections-20261005-v1/identity_resource_qualified/all_screened_roster.tsv"
COLOC = f"{ROOT}/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv"
TIER = f"{ROOT}/GWAS/finemapping/config/gwas_trait_tier.tsv"


def bh(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    o = np.argsort(p)
    ranks = np.empty(n)
    ranks[o] = np.arange(1, n + 1)
    q = p * n / ranks
    qs = np.minimum.accumulate(q[o][::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(qs, 1)
    return out


def c1():
    d = pd.read_csv(COHORT, sep="\t")
    print(d.groupby("dataset")[["n_control", "n_disease"]].first().to_string())
    d["SE"] = d.logFC / d.t
    rows = []
    for gid, x in d.groupby("gene_id_versioned"):
        y = x.logFC.to_numpy()
        v = x.SE.to_numpy() ** 2
        k = len(y)
        w = 1 / v
        fe = (w * y).sum() / w.sum()
        Q = (w * (y - fe) ** 2).sum()
        c = w.sum() - (w ** 2).sum() / w.sum()
        tau2 = max(0.0, (Q - (k - 1)) / c)
        ws = 1 / (v + tau2)
        mu = (ws * y).sum() / ws.sum()
        qhk = (ws * (y - mu) ** 2).sum() / (k - 1)
        sehk = np.sqrt(max(qhk, 1.0) / ws.sum())
        rows.append((gid, k, mu, sehk))
    r = pd.DataFrame(rows, columns=["gene", "k", "mu", "sehk"])
    print("genes", len(r), "k=5 genes", (r.k == 5).sum())
    for m in (0.25, 0.5):
        p = np.maximum(stats.t.sf((r.mu + m) / r.sehk, 4), stats.t.cdf((r.mu - m) / r.sehk, 4))
        print(f"margin {m}: RE-DL + truncated Knapp-Hartung TOST, BH<0.05 over {len(r)} genes:", int((bh(p) < 0.05).sum()))


def c2():
    g = pd.read_csv(f"{ROOT}/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/cohort_group_counts.tsv", sep="\t")
    print(g[["dataset", "new_control", "new_disease"]].to_string(index=False))
    print("units", g.new_control.sum() + g.new_disease.sum(), "control", g.new_control.sum(), "disease", g.new_disease.sum())
    d = pd.read_csv(POOLED)
    print("genes", len(d), "canonical calls", int(((d.padj < 0.05) & (d.logFC.abs() > 0.5)).sum()))
    print("SE quantiles 5/50/95%:", d.SE.quantile([0.05, 0.5, 0.95]).round(4).tolist())


def c3():
    U = set(pd.read_csv(POOLED).symbol.dropna())
    tot = 0
    totu = 0
    for k in ["govaere_2020", "hoang_2019", "pantano_2021", "suppli_2019", "moylan_2014", "arendt_2015"]:
        p = pd.read_csv(f"{ROOT}/data/published_gene_panels/{k}_panel.tsv", sep="\t", comment="#").dropna(subset=["gene_symbol"])
        n = len(p)
        nu = int(p.gene_symbol.isin(U).sum())
        tot += n
        totu += nu
        print(k, "claims", n, "in 24,115-gene universe", nu)
    print("TOTAL claims", tot, "testable", totu)
    for f in ["table_s2_inflammation_vs_normal.csv", "table_s3_fibrosis_vs_normal.csv", "table3_fibrosis_vs_inflammation.csv", "replicated_16_genes.csv"]:
        print(f, len(pd.read_csv(f"{ROOT}/data/published_degs/PRJNA512027/{f}")))


def c4():
    r = pd.read_csv(ROSTER, sep="\t", low_memory=False)
    m = pd.read_csv(f"{ROOT}/RNA-seq/Human/Patient_Cohorts/pipelines/custom/PRJNA512027/metadata/SraRunTable.csv")
    m["batch"] = m["Library Name"].str[:2]
    x = r[(r.dataset == "PRJNA512027") & (r.primary_selected == True)].merge(m[["Run", "batch"]], left_on="sample_id", right_on="Run")
    print(pd.crosstab(x.condition, x.batch).to_string())
    y = r[(r.dataset == "GSE240729") & (r.primary_selected == True)]
    print("GSE240729 stage counts", y.fibrosis_stage.value_counts().sort_index().to_dict())


def c5():
    t = pd.read_csv(TIER, sep="\t")
    t12 = set(t[t.tier.astype(str).isin(["1", "2"])].study_name)
    a = pd.read_csv(COLOC, low_memory=False, usecols=["gwas_name", "ensembl", "PP.H3.susie", "PP.H4.susie", "method"])
    s = a[(a.method == "susie") & (a.gwas_name.isin(t12))]
    g = s.groupby("ensembl").agg(h4=("PP.H4.susie", "max"), h3=("PP.H3.susie", "max"))
    print("tier-1/2 GWAS", len(t12), "SuSiE-tested genes", len(g), "max H4<=0.5", int((g.h4 <= 0.5).sum()), "and max H3>0.8", int(((g.h4 <= 0.5) & (g.h3 > 0.8)).sum()))
    s = a[a.method == "susie"].drop_duplicates(["gwas_name", "ensembl"])
    for A, B in [("UKBB_ALT", "MVP_ALT_EUR"), ("UKBB_AST", "MVP_AST_EUR")]:
        n = s[s.gwas_name == A][["ensembl"]].merge(s[s.gwas_name == B][["ensembl"]], on="ensembl")
        print(A, B, "genes with SuSiE-COLOC in both", len(n))


def c6():
    t = pd.read_csv(TIER, sep="\t")
    t12 = set(t[t.tier.astype(str).isin(["1", "2"])].study_name)
    a = pd.read_csv(COLOC, low_memory=False, usecols=["gwas_name", "ensembl", "method", "PP.H4.susie"])
    a = a[a.gwas_name.isin(t12)]
    print("tier-1/2 gene-GWAS pairs", len(a), "SuSiE pairs", int((a.method == "susie").sum()), "SuSiE PP.H4>0.9 pairs", int((a["PP.H4.susie"] > 0.9).sum()))


if __name__ == "__main__":
    {"C1": c1, "C2": c2, "C3": c3, "C4": c4, "C5": c5, "C6": c6}[sys.argv[1]]()
