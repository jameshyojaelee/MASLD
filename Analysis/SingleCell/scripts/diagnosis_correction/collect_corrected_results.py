"""Collect corrected estimates and their denominators for review text."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2

root = Path(__file__).resolve().parents[4]
task = root / "Analysis/SingleCell/candidates/source-diagnosis-correction-20261026"


def read(name):
    return pd.read_csv(task / name, sep="\t")


def record(frame):
    return json.loads(frame.to_json(orient="records"))


meta = read("metadata/donor_metadata.tsv")
p = read("programs/results_primary_family.tsv")
qc = read("programs_qc/results_primary_family.tsv")
ambient = read("programs_ambient_v2/results_primary_family.tsv")
tc = read("programs_tcell20/results_primary_family.tsv")
assert len(p) == 117 and len(ambient) == 106
v = p.copy()
v["cohort_Q"] = (v.slope_GSE202379 - v.slope_GSE244832)**2 / (v.se_GSE202379**2 + v.se_GSE244832**2)
v["cohort_Q_p"] = chi2.sf(v.cohort_Q, 1)
order = v.cohort_Q_p.sort_values().index
v.loc[order, "cohort_Q_q"] = np.minimum.accumulate((v.loc[order, "cohort_Q_p"].to_numpy() * 117 / np.arange(1, 118))[::-1])[::-1].clip(0, 1)
scores = read("programs/donor_program_scores_primary.tsv")
fractions = []
for _, row in p[p.cell_type.eq("hepatocytes") & p.module.isin([8, 20])].iterrows():
    d = scores[scores.cell_type.eq(row.cell_type) & scores.module.eq(row.module)]
    information = []
    for ds in ["GSE202379", "GSE244832"]:
        stages = d.loc[d.dataset.eq(ds), "stage"]
        information.append(((stages - stages.mean())**2).sum() / row["sigma_" + ds]**2)
    fractions.append(dict(program=row.program, GSE202379_information_fraction=information[0] / sum(information)))
g = pd.read_csv(task / "pseudobulk_de_v2/Hepatocytes_de.csv")
strict = pd.read_csv(task / "pseudobulk_de_strict_control/Hepatocytes_de.csv")
roster = read("pseudobulk_de_v2/Hepatocytes_donor_roster.tsv")
lineage = read("content_lineage/lineage_meta_effects.tsv")
lineage_primary = lineage[lineage.minimum_cells.eq(50) & lineage.annotation_filter.eq("all") &
    lineage.arm.eq("ambient_corrected_original_exposure") & lineage.contrast.eq("primary_steatohepatitis_vs_steatosis")]
ccc = read("CCC/fig4d_program_ccc_associations.tsv")
full_ccc = read("CCC/full_donor_lr_family_dispositions.tsv")
fixed_ccc = read("CCC/fixed_13_current_stage_context.tsv")
output = dict(
    status="corrected candidate; not adopted into named manuscript release",
    descriptive_diagnoses=record(meta[meta.analyzed_donor].groupby("disease_stage_coarse").size().rename("n").reset_index()),
    primary_n=record(p.groupby("cell_type").first().reset_index()[["cell_type", "n_fit", "n_Healthy", "n_Steatosis", "n_Steatohepatitis"]]),
    programs=dict(selected=int(p.selected_H59.sum()), smallest_q=float(p.perm_q.min()),
        focal=record(v[v.cell_type.eq("hepatocytes") & v.module.isin([8,20])]),
        cohort_opposite_signs=int((p.slope_GSE202379 * p.slope_GSE244832 < 0).sum()),
        cohort_Q_BH05=int(v.cohort_Q_q.lt(.05).sum()),
        cohort_smallest_Q_q=float(v.cohort_Q_q.min()), information=fractions),
    qc=dict(selected=int(qc.selected_H59.sum()), smallest_q=float(qc.perm_q.min()),
        beta_rho=float(p.set_index("program").gls_beta.corr(qc.set_index("program").gls_beta,method="spearman")),
        focal=record(qc[qc.cell_type.eq("hepatocytes") & qc.module.isin([8,20])])),
    ambient=dict(selected=int(ambient.selected_H59.sum()),smallest_q=float(ambient.perm_q.min()),
        focal=record(ambient[ambient.cell_type.eq("hepatocytes") & ambient.module.isin([8,20])])),
    tcell20=dict(n=int(tc.n_fit.iloc[0]),smallest_q=float(tc.perm_q.min()),selected=int(tc.selected_H59.sum())),
    genes=dict(hepatocyte_tested=len(g),hepatocyte_BH05=int(g.padj.lt(.05).sum()),
        hepatocyte_donors=len(roster),hepatocyte_strict_BH05=int(strict.padj.lt(.05).sum())),
    ambient_paired_focal=record(read("ambient_program_effects.tsv").query(
        "cell_type == 'hepatocytes' and module in [8,20] and analysis_universe == 'complete_case_common_universe'")),
    content_lineage=record(lineage_primary),
    communication=dict(planned=26,estimable=int(ccc.testable.sum()),supported=record(ccc[ccc.evidence_state.eq("supported")]),
        native_groups=len(full_ccc),native_fitted=int(full_ccc.disposition.eq("fitted").sum()),
        fixed_display=record(fixed_ccc)),
    core=json.loads((task / "core_verification.json").read_text()),
    models=json.loads((task / "models_verification.json").read_text()),
)
out = task / "corrected_results.json"
assert not out.exists()
out.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")
v.to_csv(task / "program_cohort_heterogeneity.tsv",sep="\t",index=False)
print("Collected corrected results and actual denominators.")
