"""F7b repair, R08-02/R08-03 and conflict 1: re-derive the Figure 7 severity-transfer values
in the primary inferential units straight from the released-scorer outputs (no refit), and
build the panel 7D table whose block b uses those units.

Primary units:
  GSE268273  109 participant records (source-declared, not genotype-screened)
  GSE213621  263 genetic identities (primary-selected record; 11 identities with two
             different stage labels removed)
  GSE276114  81 participants with MASLD; the 82 with chronic viral hepatitis are a separate row
Record-level published values are recomputed alongside as the sensitivity set.

95% intervals: percentile bootstrap, 10,000 draws, seed 20260923. Primary rows hold one
record per unit, so units are resampled within stage strata, the method of the published
intervals (external_validation.py:235-238, run_tabicl.py:206-209). The published-record
sensitivity rows of GSE213621 resample genetic-identity clusters without stratification.
"""
import json
import numpy as np
import pandas as pd
from scipy import stats

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
EX = f"{ROOT}/Analysis/MASLD_Model_Benchmark/executions"
EV = f"{EX}/external-validation-expansion-20260901T190000Z"
OUT = f"{ROOT}/figures/candidates/review-r3-fixes-20261015T112148Z/F7b/repair"
F7A_TAB = f"{ROOT}/figures/candidates/review-r3-fixes-20261015T112148Z/F7a/source_tables/panel_7d_bounds.tsv"
SEED, B = 20260923, 10000
rng = np.random.default_rng(SEED)


def tie_ceiling(y):  # verbatim from external_validation.py (pipeline copy of job2_ood_transfer.py)
    y = np.asarray(y, float)
    if len(np.unique(y)) < 2:
        return float("nan")
    mr = stats.rankdata(y, method="average")
    return float(np.corrcoef(np.sort(mr), np.arange(1, len(y) + 1))[0, 1])


def rho(s, y):
    return float(stats.spearmanr(s, y).statistic)


def unit_boot(d, unit):
    """One record per unit: stage-stratified resampling of units. Several records per unit
    (published GSE213621 record sets): unstratified resampling of identity clusters."""
    out = np.empty(B)
    sc, yy = d.score.to_numpy(float), d.y.to_numpy(float)
    if d[unit].is_unique:
        strata = [np.where(yy == u)[0] for u in np.unique(yy)]
        for b in range(B):
            take = np.concatenate([rng.choice(st, size=len(st), replace=True) for st in strata])
            out[b] = rho(sc[take], yy[take])
        how = "stage-stratified unit"
    else:
        idx = d.groupby(unit).indices
        keys = np.array(list(idx))
        for b in range(B):
            take = np.concatenate([idx[k] for k in rng.choice(keys, len(keys))])
            out[b] = rho(sc[take], yy[take])
        how = "identity cluster"
    return float(np.quantile(out, 0.025)), float(np.quantile(out, 0.975)), how


rows = []


def add(cohort, subset, unit_name, d, unit_col, primary):
    r = rho(d.score, d.y)
    c = tie_ceiling(d.y)
    lo, hi, how = unit_boot(d.reset_index(drop=True), unit_col)
    rows.append(dict(cohort=cohort, subset=subset, primary=primary, unit=unit_name,
                     n_units=int(d[unit_col].nunique()), n_records=len(d), spearman=r,
                     tie_ceiling=c, ratio=r / c, ci_lo=lo, ci_hi=hi, ci_B=B, ci_seed=SEED,
                     ci_unit=unit_col, ci_method=how))
    print(f"{cohort:10s} {subset:60s} n={len(d):4d} rho={r:.6f} ceil={c:.6f} ratio={r/c:.6f} "
          f"CI [{lo:.4f}, {hi:.4f}]", flush=True)


# ---------------------------------------------------------------- GSE268273
ph = pd.read_csv(f"{EX}/gse268273-evaluator-phenotype-20260830T205701Z/gse268273_participant_phenotype.tsv",
                 sep="\t", dtype={"row_id": str})
pr = pd.read_csv(f"{EX}/release-instrument-20260901T211017Z/sealed/jobC_gse268273_predictions.tsv",
                 sep="\t", dtype={"row_id": str})
exp = json.load(open(f"{ROOT}/Analysis/MASLD_Model_Benchmark/release/masld-severity-v1.2/expected_ood.json"))
arm = exp["arm"]
assert arm == "release_26629__quantile__latent"
a = pr[["row_id", arm]].merge(ph[["row_id", "fibrosis_stage"]], on="row_id", how="inner", validate="1:1")
assert len(a) == 109
a = a.rename(columns={arm: "score"})
a["y"] = a.fibrosis_stage.astype(int)
assert abs(rho(a.score, a.y) - exp["ood_spearman"]) < 1e-9
assert abs(tie_ceiling(a.y) - exp["tie_ceiling"]) < 1e-9
add("GSE268273", "fibrosis 0-4, all participant records", "participant record", a, "row_id", True)

# ---------------------------------------------------------------- GSE213621
s = pd.read_csv(f"{EV}/GSE213621_scores.tsv", sep="\t")
m = pd.read_csv(f"{ROOT}/Analysis/Deconvolution/bulk/GSE213621/GSE213621_metadata.tsv", sep="\t")
ro = pd.read_csv(f"{ROOT}/figures/candidates/resource-corrections-20261005-v1/identity_resource_qualified/all_screened_roster.tsv",
                 sep="\t", low_memory=False)
assert not ro.groupby("screen_component").dataset.nunique().gt(1).any(), "identity spans datasets"
ro = ro[ro.dataset == "GSE213621"][["sample_id", "screen_component", "primary_selected"]]
d = (s.merge(m[["sample_id", "fibrotic_stage"]], on="sample_id", validate="1:1")
      .merge(ro, on="sample_id", how="left", validate="1:1"))
assert d.screen_component.notna().all() and len(d) == 367
d["score"] = d.latent_severity
for lev, mp in [("3-level", {"F0F1": 0, "F2": 1, "F3F4": 2}),
                ("4-level", {"Control": 0, "F0F1": 1, "F2": 2, "F3F4": 3})]:
    x = d[d.fibrotic_stage.isin(mp)].copy()
    x["y"] = x.fibrotic_stage.map(mp)
    disc = x.groupby("screen_component").y.nunique()
    bad = set(disc[disc > 1].index)
    multi = x.screen_component.value_counts()
    print(lev, "records", len(x), "identities", x.screen_component.nunique(),
          "two-library identities", int((multi > 1).sum()), "discordant", len(bad))
    p = x[(x.primary_selected == True) & ~x.screen_component.isin(bad)]  # noqa: E712
    assert p.screen_component.is_unique
    add("GSE213621", f"{lev}, one record per genetic identity, discordant removed", "genetic identity",
        p, "screen_component", lev == "3-level")
    add("GSE213621", f"{lev}, all records (published set)", "record (identity-clustered CI)",
        x, "screen_component", False)

# ---------------------------------------------------------------- GSE276114
s = pd.read_csv(f"{EV}/GSE276114_scores.tsv", sep="\t")
dz = pd.read_csv(f"{EX}/etiology-transfer-fibrosis-distinctness-20260901T135543Z/gse276114_design.tsv", sep="\t")
d = s.merge(dz, left_on="sample_id", right_on="matrix_column_name", validate="1:1")
assert len(d) == 177
d["y"] = d.disease_group.map({"F0-2": 0, "F3": 1, "F4": 2})
assert d.y.notna().all()
d["score"] = d.latent_severity
print(pd.crosstab(d.disease, d.disease_group))
add("GSE276114", "MASLD, 3-level stage", "participant record", d[d.disease == "MASLD"], "sample_id", True)
add("GSE276114", "chronic viral hepatitis, 3-level stage", "participant record", d[d.disease == "CVH"],
    "sample_id", True)
add("GSE276114", "all etiologies, 3-level stage (published set)", "participant record", d, "sample_id", False)
add("GSE276114", "CVH plus ARLD, 3-level stage", "participant record", d[d.disease != "MASLD"],
    "sample_id", False)

t = pd.DataFrame(rows)
t.to_csv(f"{OUT}/fig7_transfer_primary_units.tsv", sep="\t", index=False, float_format="%.6f")
prim = t[t.primary]
print("primary ratio range", prim[prim.subset != "chronic viral hepatitis, 3-level stage"].ratio.min(),
      prim[prim.subset != "chronic viral hepatitis, 3-level stage"].ratio.max())
print("MASLD participant total",
      int(prim[prim.subset != "chronic viral hepatitis, 3-level stage"].n_units.sum()))

# ---------------------------------------------------------------- panel 7D table
BLK = ("b. released severity instrument in its three transfer cohorts: "
       "Spearman beside the tie-structure ceiling")
f7a = pd.read_csv(F7A_TAB, sep="\t", dtype=str, keep_default_na=False)
assert (f7a.block == BLK).sum() == 12
keep = f7a[f7a.block != BLK]
LBL = {
    ("GSE268273", True, "fibrosis 0-4, all participant records"): "GSE268273 fibrosis 0-4, n=109 participant records",
    ("GSE213621", True, "3-level, one record per genetic identity, discordant removed"):
        "GSE213621 3-level stage, n=263 genetic identities",
    ("GSE276114", True, "MASLD, 3-level stage"): "GSE276114 MASLD 3-level stage, n=81 participants",
    ("GSE276114", True, "chronic viral hepatitis, 3-level stage"):
        "GSE276114 chronic viral hepatitis 3-level stage, n=82 participants",
}
SRC = {
    "GSE268273": "Analysis/MASLD_Model_Benchmark/executions/release-instrument-20260901T211017Z/sealed/jobC_gse268273_predictions.tsv + gse268273-evaluator-phenotype-20260830T205701Z/gse268273_participant_phenotype.tsv",
    "GSE213621": "Analysis/MASLD_Model_Benchmark/executions/external-validation-expansion-20260901T190000Z/GSE213621_scores.tsv + GSE213621_metadata.tsv + identity_resource_qualified/all_screened_roster.tsv",
    "GSE276114": "Analysis/MASLD_Model_Benchmark/executions/external-validation-expansion-20260901T190000Z/GSE276114_scores.tsv + etiology-transfer-fibrosis-distinctness-20260901T135543Z/gse276114_design.tsv",
}
new = []
for _, r in prim.iterrows():
    lab = LBL[(r.cohort, True, r.subset)]
    note = f"released v1.2 weights, no refit; ratio to ceiling {r.ratio:.4f}"
    for q, v, k in [("achieved Spearman (released v1.1/v1.2 weights, no refit)", r.spearman, "achieved"),
                    ("tie-structure ceiling", r.tie_ceiling, "ceiling"),
                    ("achieved Spearman 95% CI low", r.ci_lo, "ci_lo"),
                    ("achieved Spearman 95% CI high", r.ci_hi, "ci_hi")]:
        new.append(dict(block=BLK, row=lab, quantity=q, value=repr(float(v)), kind=k, note=note,
                        source=SRC[r.cohort]))
assert len(new) == 16
out = pd.concat([keep, pd.DataFrame(new)], ignore_index=True)
# keep F7a's block order a, b, c, d
out["_o"] = out.block.str[0]
out = out.sort_values("_o", kind="stable").drop(columns="_o")
import os
os.makedirs(f"{OUT}/tables", exist_ok=True)
out.to_csv(f"{OUT}/tables/panel_7d_bounds.tsv", sep="\t", index=False)
print("wrote panel table rows", len(out))
