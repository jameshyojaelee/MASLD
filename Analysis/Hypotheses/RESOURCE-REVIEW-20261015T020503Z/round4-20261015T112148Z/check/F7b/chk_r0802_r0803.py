"""Blind check of R08-02 / R08-03 numbers from release-level files.
Independent code path: tie ceiling computed as the Spearman of y with a tie-broken
perfect ordering; identity clusters from the roster screen_component; one record per
cluster chosen two ways (roster primary_selected, and a deterministic rule that does
not read primary_selected). Seed 20260923."""
import json
import numpy as np, pandas as pd
from scipy import stats

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
E = f"{ROOT}/Analysis/MASLD_Model_Benchmark/executions/external-validation-expansion-20260901T190000Z"
OUT = f"{ROOT}/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/F7b"
SEED = 20260923


def ceiling(y):
    y = np.asarray(y, float)
    order = np.argsort(y, kind="mergesort")
    perfect = np.empty(len(y)); perfect[order] = np.arange(len(y))
    return stats.spearmanr(perfect, y).statistic


def rho(s, y):
    return stats.spearmanr(s, y).statistic


res = {}
# ---------------- GSE276114 by etiology ----------------
sc = pd.read_csv(f"{E}/GSE276114_scores.tsv", sep="\t")
de = pd.read_csv(f"{ROOT}/Analysis/MASLD_Model_Benchmark/executions/etiology-transfer-fibrosis-distinctness-20260901T135543Z/gse276114_design.tsv", sep="\t")
d = sc.merge(de, left_on="sample_id", right_on="matrix_column_name", how="inner")
assert len(d) == 177 and d.geo_accession.is_unique
d["y"] = d.disease_group.map({"F0-2": 0, "F3": 1, "F4": 2})
assert d.y.notna().all()
res["GSE276114_crosstab"] = pd.crosstab(d.disease, d.disease_group).to_dict()
for lab, sub in [("all", d), ("MASLD", d[d.disease == "MASLD"]), ("CVH", d[d.disease == "CVH"]),
                 ("ARLD", d[d.disease == "ARLD"]), ("nonMASLD", d[d.disease != "MASLD"])]:
    if sub.y.nunique() < 2:
        res[f"GSE276114_{lab}"] = dict(n=len(sub), levels=int(sub.y.nunique()))
        continue
    r, c = rho(sub.latent_severity, sub.y), ceiling(sub.y)
    res[f"GSE276114_{lab}"] = dict(n=len(sub), rho=r, ceil=c, ratio=r / c)

# ---------------- GSE213621 by identity ----------------
sc = pd.read_csv(f"{E}/GSE213621_scores.tsv", sep="\t")
me = pd.read_csv(f"{ROOT}/Analysis/Deconvolution/bulk/GSE213621/GSE213621_metadata.tsv", sep="\t")
# metadata has a ragged N+1 header; check alignment explicitly on a known row
assert me.loc[me.sample_id == "SRR21622826", "fibrotic_stage"].iloc[0] == "Control"
ro = pd.read_csv(f"{ROOT}/figures/candidates/resource-corrections-20261005-v1/identity_resource_qualified/all_screened_roster.tsv", sep="\t", low_memory=False)
cross = ro.groupby("screen_component").dataset.nunique()
res["roster_components_spanning_datasets"] = int((cross > 1).sum())
res["roster_has_GSE268273"] = int((ro.dataset == "GSE268273").sum())
res["roster_has_GSE276114"] = int((ro.dataset == "GSE276114").sum())
r6 = ro[ro.dataset == "GSE213621"][["sample_id", "screen_component", "primary_selected", "source_fibrosis_category", "genetic_person_verified", "exclusion_reason"]]
d = sc.merge(me[["sample_id", "fibrotic_stage"]], on="sample_id", how="inner").merge(r6, on="sample_id", how="inner")
assert len(d) == 367
# stage label agreement between metadata and roster
mm = {"F0F1": "F0–F1", "F2": "F2", "F3F4": "F3–F4", "Control": "Control"}
res["meta_vs_roster_stage_disagree"] = int((d.fibrotic_stage.map(mm) != d.source_fibrosis_category).sum())
rng = np.random.default_rng(SEED)
for arm, mp in [("3level", {"F0F1": 0, "F2": 1, "F3F4": 2}), ("4level", {"Control": 0, "F0F1": 1, "F2": 2, "F3F4": 3})]:
    a = d[d.fibrotic_stage.isin(mp.keys())].copy(); a["y"] = a.fibrotic_stage.map(mp)
    sizes = a.groupby("screen_component").size()
    nlev = a.groupby("screen_component").y.nunique()
    disc = set(nlev[nlev > 1].index)
    pair_types = (a[a.screen_component.isin(disc)].groupby("screen_component").fibrotic_stage
                  .apply(lambda x: "/".join(sorted(x))).value_counts().to_dict())
    r, c = rho(a.latent_severity, a.y), ceiling(a.y)
    out = dict(n=len(a), clusters=int(a.screen_component.nunique()), multi=int((sizes > 1).sum()),
               max_size=int(sizes.max()), discordant=len(disc), discordant_types=pair_types,
               published_rho=r, published_ceil=c, published_ratio=r / c)
    # (i) roster primary_selected, discordant dropped
    p = a[(a.primary_selected == True) & ~a.screen_component.isin(disc)]
    out["primary_selected_unique_per_cluster"] = bool(p.screen_component.is_unique)
    out["primary_n"] = len(p)
    out["clusters_after_drop"] = int(a[~a.screen_component.isin(disc)].screen_component.nunique())
    rp, cp = rho(p.latent_severity, p.y), ceiling(p.y)
    out.update(primary_rho=rp, primary_ceil=cp, primary_ratio=rp / cp)
    # (ii) rule that ignores primary_selected: lowest sample_id per concordant cluster
    q = a[~a.screen_component.isin(disc)].sort_values("sample_id").drop_duplicates("screen_component")
    rq, cq = rho(q.latent_severity, q.y), ceiling(q.y)
    out.update(alt_n=len(q), alt_rho=rq, alt_ratio=rq / cq)
    # (iii) cluster-mean score
    g = a[~a.screen_component.isin(disc)].groupby("screen_component").agg(s=("latent_severity", "mean"), y=("y", "first"))
    out.update(mean_rho=rho(g.s, g.y), mean_ratio=rho(g.s, g.y) / ceiling(g.y))
    # identity bootstrap on (i): one record per identity, so resample records
    s_arr, y_arr = p.latent_severity.to_numpy(), p.y.to_numpy()
    bs = np.array([rho(s_arr[i], y_arr[i]) for i in (rng.integers(0, len(p), len(p)) for _ in range(4000))])
    out["primary_boot95"] = [float(np.quantile(bs, .025)), float(np.quantile(bs, .975))]
    # identity bootstrap on the published set (resample clusters)
    grp = {k: v.to_numpy() for k, v in a.reset_index(drop=True).groupby("screen_component").groups.items()}
    keys = list(grp); A = a.reset_index(drop=True)
    sa, ya = A.latent_severity.to_numpy(), A.y.to_numpy()
    bc = []
    for _ in range(4000):
        idx = np.concatenate([grp[keys[k]] for k in rng.integers(0, len(keys), len(keys))])
        bc.append(rho(sa[idx], ya[idx]))
    out["published_cluster_boot95"] = [float(np.quantile(bc, .025)), float(np.quantile(bc, .975))]
    res[f"GSE213621_{arm}"] = out

# GSE268273 and paired deltas from release files
eo = json.load(open(f"{ROOT}/Analysis/MASLD_Model_Benchmark/release/masld-severity-v1.2/expected_ood.json"))
res["GSE268273"] = dict(rho=eo["ood_spearman"], ceil=eo["tie_ceiling"], ratio=eo["ood_spearman"] / eo["tie_ceiling"], n=eo["n_participants"],
                        point_delta=eo["ood_spearman"] - eo["native_in_cohort"])
rj = json.load(open(f"{E}/results.json"))["results"]
res["deltas"] = {k: dict(delta=v["vs_native_paired_delta"], ci=v["vs_native_paired_delta_ci95"], mde=v["vs_native_mde_80pct"]) for k, v in rj.items()}
res["sum_453"] = 109 + res["GSE213621_3level"]["primary_n"] + res["GSE276114_MASLD"]["n"]
json.dump(res, open(f"{OUT}/chk_r0802_r0803.json", "w"), indent=1, default=float)
print(json.dumps(res, indent=1, default=float))
