"""Fig 4F recomputation: cohort-standardized Ductular injury and Stromal ECM scores (outcome_z; SD 1 within cohort) averaged
over nine 20%-wide windows of the within-cohort fixed_projection percentile, centres 10..90%. Window rule [lower, upper),
last window tested both ways. Also window n, SE, and the per-window GSE126848 sample sizes."""
import numpy as np, pandas as pd
root = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
out = f"{root}/Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig4/recompute"
HAC = f"{root}/RNA-seq/results/histology_anchored_continuum"
P = pd.read_csv(f"{HAC}/candidates/hac-continuum-corrected774-20261009T171017Z/projection/participant_scores.tsv", sep="\t")
P = P[P.axis_id == "fixed_projection"][["sample_id", "dataset", "axis_percentile", "fibrosis_stage"]]
H = pd.read_csv(f"{HAC}/molecular_layers/hac-molecular-layers-corrected774-20261009T171017Z/programs/hotspot/hotspot_participant_scores.tsv.gz", sep="\t",
                usecols=["program_uid", "sample_id", "dataset", "score_raw", "outcome_z", "module_name"])
foc = {"hotspot_hepatocytes_f05c535ae5bbc0b9": None, "hotspot_hepatocytes_48f39dd4d817a10e": None}
H = H[H.program_uid.isin(foc)]
print(H[["program_uid", "module_name"]].drop_duplicates().to_string())
D = H.merge(P, on=["sample_id", "dataset"], how="inner")
print("units per cohort", D.groupby(["program_uid", "dataset"]).size().unstack())
rows = []
for (pu, ds), g in D.groupby(["program_uid", "dataset"]):
    z = (g.score_raw - g.score_raw.mean()) / g.score_raw.std(ddof=1)
    for k, c in enumerate(np.round(np.arange(0.1, 0.91, 0.1), 2)):
        lo, hi = c - 0.1, c + 0.1
        m1 = (g.axis_percentile >= lo - 1e-12) & (g.axis_percentile < hi - 1e-12)
        m2 = m1 | ((g.axis_percentile >= lo - 1e-12) & (abs(hi - 1.0) < 1e-9) & (g.axis_percentile <= 1.0))
        for rule, m in [("half_open", m1), ("last_closed", m2)]:
            rows.append(dict(program_uid=pu, module_name=g.module_name.iloc[0], dataset=ds, centre=c, rule=rule, n=int(m.sum()),
                             mean_outcome_z=g.outcome_z[m].mean(), mean_rez=z[m].mean(), se=g.outcome_z[m].std(ddof=1) / np.sqrt(m.sum())))
W = pd.DataFrame(rows)
W.to_csv(f"{out}/r4f_windows.tsv", sep="\t", index=False)
A = pd.read_csv(f"{HAC}/molecular_layers/hac-molecular-layers-corrected774-20261009T171017Z/programs/hotspot/hotspot_fixed_windows.tsv.gz", sep="\t")
print(A.columns.tolist())
print(A.head(2).to_string())
