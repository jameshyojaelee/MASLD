"""DIG1 job 12 (exploratory, descriptive). Do receptors carried by GNMT-low hepatocytes have
ligands that rise with MASLD in another lineage, in both cohorts?
State receptor: liana consensus receptor subunit that is a negative within-donor GNMT partner
(job 9, BH < 0.05) and up in MASLD hepatocytes (GSE298719 FDR < 0.10, GSE256398 logFC > 0 and
P < 0.05). Comparison: all other liana receptor subunits tested in job 9 (detected in >= 10% of
hepatocytes). Ligand rise: logFC > 0 and P < 0.05 in the sender lineage in every cohort where that
lineage was tested (fibroblasts and cholangiocytes: both cohorts; macrophages, T cells: GSE298719
only). Fisher exact test on the two receptor groups is descriptive (genes are not independent).
Usage: python -I 12_niche_ligands.py
"""
from pathlib import Path

import liana as li
import pandas as pd
from scipy import stats

RUN = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/DIG1-20261024T113951Z")
OUT = RUN / "out_niche_ligands"
assert not any(OUT.iterdir()), "refusing to overwrite"
lr = li.resource.select_resource("consensus")
lr = lr.assign(rec=lr.receptor.str.split("_")).explode("rec").assign(lig=lambda d: d.ligand.str.split("_")).explode("lig")
part = pd.read_csv(RUN / "out_hep_axes/gnmt_within_donor_partners.tsv.gz", sep="\t").set_index("gene")
de = pd.read_csv(RUN / "out_pseudobulk_v2/gene_DE_edgeR_all_lineages.tsv.gz", sep="\t").drop_duplicates(["lkey", "symbol"])
D = {k: g.set_index("symbol") for k, g in de.groupby("lkey")}
hF, hV = D["GSE298719:hepatocytes"], D["GSE256398:hepatocytes"]
recs = sorted(set(lr.rec) & set(part.index))
r = pd.DataFrame(index=recs)
r["gnmt_z"] = part.mean_z.reindex(recs); r["gnmt_q"] = part.bh_q.reindex(recs)
r["lfc_F"] = hF.logFC.reindex(recs); r["fdr_F"] = hF.FDR.reindex(recs)
r["lfc_V"] = hV.logFC.reindex(recs); r["p_V"] = hV.P.reindex(recs)
r["state"] = (r.gnmt_q < 0.05) & (r.gnmt_z < 0) & (r.lfc_F > 0) & (r.fdr_F < 0.10) & (r.lfc_V > 0) & (r.p_V < 0.05)
SENDERS = {"fibroblasts": ["GSE298719:fibroblasts", "GSE256398:fibroblasts"],
           "cholangiocytes": ["GSE298719:cholangiocytes", "GSE256398:cholangiocytes"],
           "macrophages": ["GSE298719:macrophages"], "tcells": ["GSE298719:tcells"],
           "hepatocytes": ["GSE298719:hepatocytes", "GSE256398:hepatocytes"]}
rows = []
for rec in recs:
    for lig in sorted(set(lr.lig[lr.rec == rec])):
        for s, keys in SENDERS.items():
            vals = [(D[k].logFC.get(lig), D[k].P.get(lig)) for k in keys]
            if any(v[0] is None or pd.isna(v[0]) for v in vals):
                continue
            rows.append(dict(receptor=rec, ligand=lig, sender=s, n_cohorts=len(keys),
                             rise=all(l > 0 and p < 0.05 for l, p in vals),
                             logFC=";".join(f"{l:.2f}" for l, _ in vals), P=";".join(f"{p:.1e}" for _, p in vals)))
L = pd.DataFrame(rows)
L = L.merge(r[["state"]], left_on="receptor", right_index=True)
L.to_csv(OUT / "receptor_ligand_sender_all.tsv", sep="\t", index=False)
any_rise = L[L.sender != "hepatocytes"].groupby("receptor").rise.any()
any_rise2 = L[(L.sender != "hepatocytes") & (L.n_cohorts == 2)].groupby("receptor").rise.any()
r["paracrine_rise"] = any_rise.reindex(recs).fillna(False)
r["paracrine_rise_2cohort"] = any_rise2.reindex(recs).fillna(False)
r.to_csv(OUT / "hepatocyte_receptors.tsv", sep="\t")
summ = []
for col in ["paracrine_rise", "paracrine_rise_2cohort"]:
    t = pd.crosstab(r.state, r[col]).reindex(index=[True, False], columns=[True, False], fill_value=0)
    orr, p = stats.fisher_exact(t.to_numpy())
    summ.append(dict(test=col, state_with=t.loc[True, True], state_n=t.loc[True].sum(),
                     other_with=t.loc[False, True], other_n=t.loc[False].sum(), odds_ratio=orr, fisher_p=p))
pd.DataFrame(summ).to_csv(OUT / "enrichment.tsv", sep="\t", index=False)
print(pd.DataFrame(summ).to_string())
print(r[r.state].sort_values("gnmt_z").round(3).to_string())
print(L[L.state & L.rise & (L.sender != "hepatocytes")].sort_values(["receptor", "sender"]).to_string())
