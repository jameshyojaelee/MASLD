"""F7b blind re-check: re-derive R08-02/03/07/08/09/10 and R11-01 numbers from producing files.
Own code, seed 20260923. Writes into the directory given as argv[1]."""
import json
import sys

import numpy as np
import pandas as pd
from scipy import stats

OUT = sys.argv[1]
ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
EX = f"{ROOT}/Analysis/MASLD_Model_Benchmark/executions"
SEED = 20260923
B = 10000
res = {}


def ceiling(y):
    """Largest Spearman any untied score can reach against y (ties in y kept)."""
    r = stats.rankdata(y)
    return float(np.corrcoef(np.sort(r), np.arange(1, len(y) + 1))[0, 1])


def strat_boot(score, y, unit_ids, rng):
    """Resample units with replacement within each stage level; one row per unit here."""
    score = np.asarray(score, float); y = np.asarray(y, float)
    groups = [np.flatnonzero(y == lv) for lv in np.unique(y)]
    out = np.empty(B)
    for b in range(B):
        idx = np.concatenate([g[rng.integers(0, len(g), len(g))] for g in groups])
        out[b] = stats.spearmanr(score[idx], y[idx]).statistic
    return np.quantile(out, [0.025, 0.975])


def summarise(name, score, y, units, rng):
    rho = stats.spearmanr(score, y).statistic
    c = ceiling(y)
    lo, hi = strat_boot(score, y, units, rng)
    r = dict(n=int(len(y)), n_unique_units=int(pd.Series(units).nunique()), spearman=float(rho),
             tie_ceiling=c, ratio=float(rho / c), ci_lo=float(lo), ci_hi=float(hi))
    res[name] = r
    print(name, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()}, flush=True)


rng = np.random.default_rng(SEED)

# ---------- GSE268273 ----------
pr_ = pd.read_csv(f"{EX}/release-instrument-20260901T211017Z/sealed/jobC_gse268273_predictions.tsv", sep="\t")
ph = pd.read_csv(f"{EX}/gse268273-evaluator-phenotype-20260830T205701Z/gse268273_participant_phenotype.tsv", sep="\t")
d = pr_.merge(ph[["row_id", "patient_code", "fibrosis_stage"]], on="row_id", validate="1:1")
assert len(d) == 109
print("GSE268273 patients", d.patient_code.nunique(), "stage counts", d.fibrosis_stage.value_counts().sort_index().to_dict())
summarise("GSE268273", d["release_26629__quantile__latent"], d.fibrosis_stage, d.patient_code, rng)
arms = ["release_26629__prospective__latent", "release_26629__prospective__ridge",
        "release_26629__quantile__latent", "release_26629__quantile__ridge"]
allarms = {c: float(stats.spearmanr(d[c], d.fibrosis_stage).statistic) for c in pr_.columns if c != "row_id"}
res["GSE268273_arm_spearman"] = allarms
print("eligible arms", {a: round(allarms[a], 4) for a in arms})
print("all arms", {a: round(v, 4) for a, v in allarms.items()})

# ---------- GSE213621 ----------
s = pd.read_csv(f"{EX}/external-validation-expansion-20260901T190000Z/GSE213621_scores.tsv", sep="\t")
m = pd.read_csv(f"{ROOT}/Analysis/Deconvolution/bulk/GSE213621/GSE213621_metadata.tsv", sep="\t")
ro = pd.read_csv(f"{ROOT}/figures/candidates/resource-corrections-20261005-v1/identity_resource_qualified/all_screened_roster.tsv",
                 sep="\t", low_memory=False)
ro = ro[ro.dataset == "GSE213621"][["sample_id", "screen_component", "primary_selected"]]
g = s.merge(m[["sample_id", "fibrotic_stage"]], on="sample_id", validate="1:1").merge(ro, on="sample_id", how="left", validate="1:1")
print("GSE213621 scored", len(g), "missing identity", int(g.screen_component.isna().sum()))
mp3 = {"F0F1": 0, "F2": 1, "F3F4": 2}
a = g[g.fibrotic_stage.isin(mp3)].copy(); a["y"] = a.fibrotic_stage.map(mp3)
nst = a.groupby("screen_component").y.nunique()
disc = set(nst[nst > 1].index)
print("3-level records", len(a), "identities", a.screen_component.nunique(), "discordant identities", len(disc),
      "multi-record identities", int((a.groupby('screen_component').size() > 1).sum()))
keep = a[~a.screen_component.isin(disc)]
# one record per identity: the roster's primary-selected record
p = keep[keep.primary_selected.astype(str).str.lower().isin(["true", "1"])]
print("primary-selected after discordant removal", len(p), "unique identities", p.screen_component.nunique(),
      "identities kept", keep.screen_component.nunique())
summarise("GSE213621_identity", p.latent_severity, p.y, p.screen_component, rng)
res["GSE213621_counts"] = dict(records_3level=int(len(a)), identities_3level=int(a.screen_component.nunique()),
                               discordant=len(disc), used=int(len(p)))
# published all-records set point estimate for reference
res["GSE213621_all_records_rho"] = float(stats.spearmanr(a.latent_severity, a.y).statistic)

# ---------- GSE276114 ----------
s = pd.read_csv(f"{EX}/external-validation-expansion-20260901T190000Z/GSE276114_scores.tsv", sep="\t")
dz = pd.read_csv(f"{EX}/etiology-transfer-fibrosis-distinctness-20260901T135543Z/gse276114_design.tsv", sep="\t")
d = s.merge(dz, left_on="sample_id", right_on="matrix_column_name", validate="1:1")
d["y"] = d.disease_group.map({"F0-2": 0, "F3": 1, "F4": 2})
print(pd.crosstab(d.disease, d.disease_group))
for dis in ["MASLD", "CVH"]:
    sub = d[d.disease == dis]
    summarise(f"GSE276114_{dis}", sub.latent_severity, sub.y, sub.sample_id, rng)

# ---------- R08-07 ----------
rj = json.load(open(f"{EX}/external-validation-expansion-20260901T190000Z/results.json"))
for k in ["GSE213621__primary_3level_no_control", "GSE276114__primary_3level"]:
    r = rj["results"][k]
    res[f"R0807_{k}"] = dict(delta=r["vs_native_paired_delta"], ci=r["vs_native_paired_delta_ci95"], mde=r["vs_native_mde_80pct"])
res["R0807_ref_point"] = rj["reference_cohort_gse268273"]["transfer_spearman"] - rj["reference_cohort_gse268273"]["native_in_cohort"]
print("R08-07", {k: v for k, v in res.items() if k.startswith("R0807")})

# ---------- R08-08 ----------
sm = pd.read_csv(f"{EX}/encoder-field-score-21314196/evaluation/study_metrics.tsv", sep="\t")
pt = sm.set_index(["model_id", "held_out_study"]).donor_class_balanced_macro_f1
clean = ["GSE136103", "GSE174748", "GSE189600", "GSE202379", "GSE244832"]
dl = {st: pt["transcriptformer_tf_sapiens::linear", st] - pt["scvi_liver_latent::linear", st] for st in clean}
res["R0808"] = dict(per_study=dl, equal_mean=float(np.mean(list(dl.values()))), GSE189600_contrib=dl["GSE189600"] / 5)
print("R08-08", res["R0808"])

# ---------- R08-09 ----------
gr = pd.read_csv(f"{EX}/fig7-source-exposure-correction-20261014/tables/panel_7f_exposure_grid.tsv", sep="\t")
o = gr[gr.held_out_study != "GSE189600"]
res["R0809"] = dict(n_blocks=int(gr.model_id.nunique()), other_min=float(o.per_study_delta.min()),
                    other_max=float(o.per_study_delta.max()), n_pos=int((o.per_study_delta > 0).sum()), n_other=int(len(o)),
                    n_studies=int(gr.held_out_study.nunique()))
print("R08-09", res["R0809"])

# ---------- R11-01 observed ----------
ASP = ["steatosis", "ballooning", "lobular_inflammation"]
j = pd.read_csv(f"{EX}/gse267145-authoritative-join-21064930/participant_join.tsv", sep="\t")
L = j[ASP].astype(float).to_numpy(); n = len(L); assert n == 99 and j.participant_id.nunique() == 99


def prt(Mx):
    ev = np.linalg.eigvalsh(Mx); return float(ev.sum() ** 2 / (ev ** 2).sum())


Rl = stats.spearmanr(L).statistic
rv = pd.read_csv(f"{EX}/model-check-1044-21181334-aspect-separability/separability/gene_association_vectors.tsv", sep="\t")
Rv = stats.spearmanr(rv[ASP].to_numpy()).statistic
iu = np.triu_indices(3, 1)
res["R1101_obs"] = dict(label_rho=Rl[iu].tolist(), label_pr=prt(Rl), rna_c3=Rv[iu].tolist(), rna_pr=prt(Rv), n_genes=len(rv))
print("R11-01 observed", res["R1101_obs"])

# ---------- R11-01 independent simulation ----------
# Own construction: features x donors; per-feature association = Pearson of feature ranks with rank-averaged grade ranks.
Lr = np.column_stack([stats.rankdata(L[:, k]) for k in range(3)])
Lz = (Lr - Lr.mean(0)) / Lr.std(0)
shared = stats.rankdata(L.sum(1)); shared = (shared - shared.mean()) / shared.std()
resid = np.empty_like(Lz)
for k in range(3):
    others = np.delete(Lz, k, axis=1)
    X = np.column_stack([np.ones(n), others])
    beta = np.linalg.lstsq(X, Lz[:, k], rcond=None)[0]
    e = Lz[:, k] - X @ beta
    resid[:, k] = (e - e.mean()) / e.std()


def assoc_stats(E):
    R = stats.rankdata(E, axis=1)
    R = (R - R.mean(1, keepdims=True)) / R.std(1, keepdims=True)
    A = R @ Lz / n
    C = stats.spearmanr(A).statistic
    return float(np.abs(C[iu]).min()), prt(C)


srng = np.random.default_rng(SEED + 1)
rows = []
for G in [30605, 96460]:
    for model in ["noise", "shared", "grade_following", "grade_specific", "mixture"]:
        grid = [(0.0, 0.0)] if model == "noise" else [(f, b) for f in [0.02, 0.05, 0.1, 0.2, 0.4] for b in [0.5, 1.0]]
        for f, b in grid:
            for rep in range(2):
                E = srng.standard_normal((G, n)).astype(np.float64)
                k = int(round(f * G))
                if k:
                    sgn = srng.choice([-1.0, 1.0], k)[:, None]
                    asp = np.arange(k) % 3
                    if model == "shared":
                        D = np.broadcast_to(shared, (k, n))
                    elif model == "grade_following":
                        D = Lz[:, asp].T
                    elif model == "grade_specific":
                        D = resid[:, asp].T
                    else:
                        h = k // 2
                        D = np.vstack([np.broadcast_to(shared, (h, n)), resid[:, asp[h:]].T])
                    E[:k] += b * sgn * D
                c3, p = assoc_stats(E)
                rows.append(dict(G=G, model=model, frac=f, beta=b, rep=rep, c3=c3, pr=p))
        print("sim", G, model, "done", flush=True)
sim = pd.DataFrame(rows)
sim.to_csv(f"{OUT}/recheck_r1101_sim_draws.tsv", sep="\t", index=False, float_format="%.6f")
rng_t = sim.groupby(["G", "model"]).agg(c3_lo=("c3", "min"), c3_hi=("c3", "max"), pr_lo=("pr", "min"), pr_hi=("pr", "max")).reset_index()
f20 = sim[sim.frac == 0.2].groupby(["G", "model"]).agg(c3_lo=("c3", "min"), c3_hi=("c3", "max"), pr_lo=("pr", "min"), pr_hi=("pr", "max")).reset_index()
rng_t.to_csv(f"{OUT}/recheck_r1101_sim_ranges.tsv", sep="\t", index=False, float_format="%.4f")
f20.to_csv(f"{OUT}/recheck_r1101_sim_frac020.tsv", sep="\t", index=False, float_format="%.4f")
pd.set_option("display.width", 200)
print(rng_t.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
print(f20.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
json.dump(res, open(f"{OUT}/recheck_values.json", "w"), indent=1, default=float)
print("done")
