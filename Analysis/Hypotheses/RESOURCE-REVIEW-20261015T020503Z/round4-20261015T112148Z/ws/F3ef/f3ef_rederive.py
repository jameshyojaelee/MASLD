# F3ef re-derivation: Fig 3E/3F sensitivity values from saved refits (no new model fits).
# Reads: released 3F source table, verify2 job 4211 (independent adjacent refits), reviewer job 1721
# (z-within-cohort complexity refits), job 4191 (F1 vs F0 leave-one-cohort-out), f0-arms r8 (disease-only F0).
import pandas as pd, numpy as np, glob, os, json
ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
RUN = ROOT + "/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z"
OUT = ROOT + "/figures/candidates/review-r3-fixes-20261015T112148Z/F3ef"
WS = RUN + "/round4-20261015T112148Z/ws/F3ef"
TR = ["F0_to_F1", "F1_to_F2", "F2_to_F3", "F3_to_F4"]

rel3f = pd.read_csv(ROOT + "/figures/main/fig3_bulk_transcriptomics/source_tables/fig3f_stage_genetics_matrix.tsv", sep="\t")
adj = pd.read_csv(ROOT + "/figures/candidates/bulk-release-corrected-controls-20261008T162308Z/analysis/stage_extensions/stage_all_gene_results.tsv",
                  sep="\t", usecols=["gene_id_versioned", "gene_id_base", "gene_name", "logFC", "FDR", "contrast"])
adj = adj[adj.contrast.isin(TR)]
v2b = dict(zip(adj.gene_id_versioned, adj.gene_id_base)); v2n = dict(zip(adj.gene_id_versioned, adj.gene_name))
coloc = pd.read_csv(ROOT + "/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv")
pos = set(coloc.loc[coloc.coloc_best_susie_pp4 > 0.5, "ensembl"]); assert len(pos) == 436

# released 3F reproduces?
rtabs = {c: adj[adj.contrast == c][["gene_id_versioned", "logFC", "FDR"]] for c in TR}
def select(tabs):
    sel = set()
    for c in TR:
        t = tabs[c]; b = t.gene_id_versioned.map(v2b)
        sel |= set(b[(t.FDR < 0.05) & b.isin(pos)])
    marks = {c: int(((tabs[c].FDR < 0.05) & tabs[c].gene_id_versioned.map(v2b).isin(sel)).sum()) for c in TR}
    return sel, marks
rsel, rmarks = select(rtabs)
assert len(rsel) == 54 and sum(rmarks.values()) == 59, (len(rsel), rmarks)
assert set(rel3f.gene_id_base) == rsel

# sensitivity fits per adjacent contrast
J4211 = RUN + "/jobs/4211/fits"
def rd(p):
    g = pd.read_csv(p, sep="\t"); return g[["gene_id_versioned", "logFC", "FDR"]]
sens = {}
for f in glob.glob(J4211 + "/*.tsv.gz"):
    lab, con = os.path.basename(f).replace(".tsv.gz", "").split("__")
    sens[("v4211_" + lab, con)] = rd(f)
J1721 = RUN + "/jobs/1721/fits"
for f in glob.glob(J1721 + "/*.genes.tsv.gz"):
    lab, con = os.path.basename(f).replace(".genes.tsv.gz", "").split("__")
    con = {"F1_vs_F0": "F0_to_F1"}.get(con, con)
    if con in TR: sens[("r1721_" + lab, con)] = rd(f)
arms = pd.read_csv(ROOT + "/figures/candidates/f0-reference-arms-r8-20261008T193309Z/analysis/f0_arms/f0_arm_results.tsv.gz", sep="\t",
                   usecols=["gene_id_versioned", "logFC", "FDR", "arm", "contrast"])
for a in ["A_all_F0", "B_disease_only_F0", "C_all_F0_plus_term"]:
    for c in TR:
        x = arms[(arms.arm == a) & (arms.contrast == c)]
        if len(x): sens[("arm_" + a, c)] = x[["gene_id_versioned", "logFC", "FDR"]]

# base reproduction checks
for k in [("v4211_base", "F0_to_F1"), ("arm_A_all_F0", "F0_to_F1")]:
    m = sens[k].merge(rtabs["F0_to_F1"], on="gene_id_versioned")
    print("reproduce", k, "n_sig", (sens[k].FDR < 0.05).sum(), "max|dlogFC|", float((m.logFC_x - m.logFC_y).abs().max()))

# per-mark table for the 54 released genes x 4 contrasts
rows = []
for c in TR:
    t = rtabs[c].copy(); t["base"] = t.gene_id_versioned.map(v2b); t = t[t.base.isin(rsel)]
    for _, r in t.iterrows():
        d = dict(gene=v2n[r.gene_id_versioned], gene_id_versioned=r.gene_id_versioned, contrast=c,
                 released_logFC=r.logFC, released_FDR=r.FDR, released_mark=r.FDR < 0.05)
        for (lab, con), g in sens.items():
            if con != c: continue
            q = g.loc[g.gene_id_versioned == r.gene_id_versioned, "FDR"]
            d["q_" + lab] = float(q.iloc[0]) if len(q) else np.nan
        rows.append(d)
pm = pd.DataFrame(rows)
pm.to_csv(OUT + "/fig3f_mark_sensitivity.tsv", sep="\t", index=False)

mk = pm[pm.released_mark]
print("released marks by contrast", mk.groupby("contrast").size().to_dict())
f01 = mk[mk.contrast == "F0_to_F1"]
for col in [c for c in pm.columns if c.startswith("q_")]:
    sub = mk[mk[col].notna()]
    if len(sub): print(f"{col}: marks kept {int((sub[col] < 0.05).sum())}/{len(sub)} by contrast",
                       sub.assign(k=sub[col] < 0.05).groupby("contrast").k.agg(['sum', 'size']).to_dict('index'))
only01 = [g for g, x in mk.groupby("gene") if set(x.contrast) == {"F0_to_F1"}]
print("genes whose only mark is F0_to_F1:", len(only01))
lostB = f01[(f01.gene.isin(only01)) & (f01["q_arm_B_disease_only_F0"] >= 0.05)]
print("only-F0F1 genes losing call under disease-only F0:", len(lostB), sorted(lostB.gene), lostB["q_arm_B_disease_only_F0"].round(4).tolist())

# reselection under joint covariates (both implementations) and under disease-only F0 arm
res = []
for lab in ["v4211_joint", "r1721_joint", "arm_B_disease_only_F0"]:
    if not all((lab, c) in sens for c in TR): print(lab, "incomplete"); continue
    sel, marks = select({c: sens[(lab, c)] for c in TR})
    b2n = {v2b[v]: v2n[v] for v in v2n}
    res.append(dict(fit=lab, genes=len(sel), effects=4 * len(sel), marks=sum(marks.values()),
                    marks_by_contrast="/".join(str(marks[c]) for c in TR), kept_of_54=len(sel & rsel),
                    left=";".join(sorted(b2n[x] for x in rsel - sel)), entered=";".join(sorted(b2n[x] for x in sel - rsel))))
rs = pd.DataFrame(res); rs.to_csv(OUT + "/fig3f_reselection_sensitivity.tsv", sep="\t", index=False)
print(rs.to_string(index=False))

# adjacent contrast totals per fit
tot = []
for (lab, con), g in sorted(sens.items()):
    tot.append(dict(fit=lab, contrast=con, n_bh05=int((g.FDR < 0.05).sum())))
for c in TR: tot.append(dict(fit="released", contrast=c, n_bh05=int((rtabs[c].FDR < 0.05).sum())))
tot = pd.DataFrame(tot); tot.to_csv(OUT + "/adjacent_contrast_counts_by_fit.tsv", sep="\t", index=False)
print(tot.pivot(index="fit", columns="contrast", values="n_bh05").to_string())

# hero genes (G2-04)
hero = pm[pm.gene.isin(["GNMT", "MAT1A", "CYP2C19"])]
print(hero.drop(columns=["gene_id_versioned"]).round(4).to_string(index=False))
