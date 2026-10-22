"""F3ef repair: summarize job 7550 cohort-omission refits (saved tables only).
Writes hero-gene q values, 3F mark checks without each cohort, and F1/F4 retention."""
import pandas as pd, numpy as np, json
ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
R = ROOT + "/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/repair/F3ef/refit-20261015T171823Z"
FIG = ROOT + "/figures/candidates/review-r3-fixes-20261015T112148Z/F3ef/repair"
rel = pd.read_csv(ROOT + "/figures/candidates/bulk-release-corrected-controls-20261008T162308Z/analysis/stage_extensions/stage_extension_all_gene_results.tsv",
                  sep="\t", usecols=["gene_id_versioned", "gene_name", "logFC", "FDR", "axis", "contrast"])
rel = rel[rel.axis.isin(["fibrosis", "fibrosis_adjacent"])]
name = dict(zip(rel.gene_id_versioned, rel.gene_name))
fc = pd.read_csv(R + "/fit_counts.tsv", sep="\t")
rp = pd.read_csv(R + "/reproduction_vs_release.tsv", sep="\t")
print(rp.to_string(index=False)); print(fc.to_string(index=False))
assert (rp.max_abs_dlogFC < 1e-8).all() and (rp.n_bh05_refit == rp.n_bh05_release).all()

def fit(f, c):
    return pd.read_csv(f"{R}/fits/{f}__{c}.tsv.gz", sep="\t")

out = {}
# F1 vs F0: released calls retained without each cohort; single-cohort sign agreement
base01 = fit("base", "F0_to_F1"); calls01 = base01[base01.FDR < 0.05]
assert len(calls01) == 1338
for cc in ["GSE130970", "GSE135251", "GSE162694"]:
    d = fit("drop_" + cc, "F0_to_F1").set_index("gene_id_versioned").loc[calls01.gene_id_versioned]
    out[f"F1_drop_{cc}_retained_of_1338"] = int(((d.FDR < 0.05) & (np.sign(d.logFC.values) == np.sign(calls01.logFC.values))).sum())
signs = []
for cc in ["GSE130970", "GSE135251", "GSE162694"]:
    a = fit("alone_" + cc, "F0_to_F1").set_index("gene_id_versioned").loc[calls01.gene_id_versioned]
    signs.append(np.sign(a.logFC.values) == np.sign(calls01.logFC.values))
out["F1_same_sign_all_three_single_cohort_fits"] = float(np.mean(np.all(signs, axis=0)))
# F4 vs F0 without GSE130970
b4 = fit("base", "F4_vs_F0"); c4 = b4[b4.FDR < 0.05]; assert len(c4) == 7043
d4 = fit("drop_GSE130970", "F4_vs_F0").set_index("gene_id_versioned")
out["F4_drop_GSE130970_total"] = int((d4.FDR < 0.05).sum())
x = d4.loc[c4.gene_id_versioned]
out["F4_drop_GSE130970_retained_of_7043"] = int(((x.FDR < 0.05) & (np.sign(x.logFC.values) == np.sign(c4.logFC.values))).sum())

# Hero genes across fits
rows = []
for g in ["GNMT", "MAT1A", "CYP2C19"]:
    gid = [k for k, v in name.items() if v == g][0]
    for c in ["F0_to_F1", "F1_to_F2", "F2_to_F3", "F3_to_F4", "F2_vs_F0", "F3_vs_F0", "F4_vs_F0"]:
        for f in ["base", "drop_GSE130970", "drop_GSE135251", "drop_GSE162694"]:
            t = fit(f, c); r = t[t.gene_id_versioned == gid].iloc[0]
            rows.append(dict(gene=g, contrast=c, fit=f, logFC=r.logFC, FDR=r.FDR))
hero = pd.DataFrame(rows)
hero.to_csv(FIG + "/hero_genes_cohort_omission.tsv", sep="\t", index=False)
print(hero.pivot_table(index=["gene", "contrast"], columns="fit", values="FDR").to_string(float_format=lambda v: f"{v:.3g}"))
late = hero[hero.contrast.isin(["F2_vs_F0", "F3_vs_F0", "F4_vs_F0"])]
for f in ["base", "drop_GSE130970", "drop_GSE135251", "drop_GSE162694"]:
    s = late[late.fit == f]
    out[f"hero_F2F4_vsF0_max_q_{f}"] = float(s.FDR.max())
    out[f"hero_F2F4_vsF0_all_negative_{f}"] = bool((s.logFC < 0).all())

# Fig 3F marks under cohort omission
m3f = pd.read_csv(ROOT + "/figures/main/fig3_bulk_transcriptomics/source_tables/fig3f_stage_genetics_matrix.tsv", sep="\t")
marks = m3f[m3f.FDR < 0.05][["gene_id_versioned", "gene_name", "contrast", "logFC", "FDR"]]
assert len(marks) == 59
mrows = []
for _, r in marks.iterrows():
    d = dict(gene=r.gene_name, contrast=r.contrast, released_q=r.FDR)
    for cc in ["GSE130970", "GSE135251", "GSE162694"]:
        t = fit("drop_" + cc, r.contrast); q = t[t.gene_id_versioned == r.gene_id_versioned].iloc[0]
        d["q_drop_" + cc] = q.FDR; d["same_sign_drop_" + cc] = bool(np.sign(q.logFC) == np.sign(r.logFC))
    mrows.append(d)
mk = pd.DataFrame(mrows)
mk.to_csv(FIG + "/fig3f_marks_cohort_omission.tsv", sep="\t", index=False)
for cc in ["GSE130970", "GSE135251", "GSE162694"]:
    k = (mk["q_drop_" + cc] < 0.05) & mk["same_sign_drop_" + cc]
    out[f"3F_marks_kept_drop_{cc}"] = mk.assign(k=k).groupby("contrast").k.sum().astype(int).to_dict()
f01 = mk[mk.contrast == "F0_to_F1"]
out["3F_F0F1_marks_min_q_drop_GSE130970"] = float(f01.q_drop_GSE130970.min())
print(json.dumps(out, indent=1))
with open(R + "/summary.json", "w") as fh:
    json.dump(out, fh, indent=1)
