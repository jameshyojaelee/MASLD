"""Light U0 check (seconds, small TSVs): Methods clustering/selection rule applied
with my own code to the producer's PLINK KING tables and callable-site counts,
with labels rebuilt from unified metadata + GSE130970 source table.
The full check (u0_roster_recount.py) recomputes KING from the VCF.
"""
import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd

R = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=True)
Q = R / "figures/candidates/resource-corrections-20261005-v1/identity_resource_qualified"
A1 = R / "Analysis/MASLD_Model_Benchmark/executions/in-biopsy-allelic-a1-full-20260924T025858Z/identity"
log = open(OUT / "u0_light_log.txt", "w")


def say(*a):
    m = " ".join(map(str, a)); print(m); log.write(m + "\n")


calls = pd.read_csv(Q / "filtered_callable_site_counts.tsv", sep="\t").set_index("sample_id")
samples = calls.index.tolist()
tabs = {lab: pd.read_csv(Q / f"{lab}_filtered_related_pairs.tsv", sep="\t") for lab in ["full", "half0", "half1"]}
key = lambda a, b: tuple(sorted((a, b)))
T = {lab: {key(a, b): (n, k) for a, b, n, k in zip(t["#IID1"], t.IID2, t.NSNP, t.KINSHIP)} for lab, t in tabs.items()}
callable_ = calls.full_called >= 1000
cand = {p for p, (n, k) in T["full"].items() if k >= 0.354 and n >= 500 and callable_[p[0]] and callable_[p[1]]}
a1 = pd.read_csv(A1 / "kinship_related_pairs.tsv", sep="\t")
a1s = {key(a, b) for a, b, r in zip(a1["#IID1"], a1.IID2, a1.relation) if r == "same_individual"}
say("candidates KING", len(cand), "A1", len(a1s), "union", len(cand | a1s), "A1-only", len(a1s - cand))
cand |= a1s


def ok(p):
    for lab, mn in [("full", 500), ("half0", 250), ("half1", 250)]:
        n, k = T[lab].get(p, (0, -1))
        if n < mn or k < 0.354:
            return False
    return bool(callable_.get(p[0], False) and callable_.get(p[1], False))


stable = {p for p in cand if ok(p)}; failed = cand - stable
say("stable", len(stable), "failed", len(failed))


def comps(nodes, edges):
    par = {s: s for s in nodes}
    def f(s):
        while par[s] != s:
            par[s] = par[par[s]]; s = par[s]
        return s
    for a, b in edges:
        ra, rb = f(a), f(b)
        if ra != rb:
            par[max(ra, rb)] = min(ra, rb)
    return {s: f(s) for s in nodes}


cs, ca = comps(samples, stable), comps(samples, stable | failed)
meta = pd.read_csv(R / "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv", keep_default_na=False, na_values=[""])
meta = meta.drop_duplicates("sample_id").set_index("sample_id")
s130 = pd.read_csv(R / "RNA-seq/Human/Patient_Cohorts/metadata/source/GSE130970/GSE130970_source_diagnosis.tsv", sep="\t").set_index("run")
L = pd.DataFrame(index=samples)
for c in ["dataset", "group_binary", "condition", "diagnosis_harmonized"]:
    L[c] = meta[c].reindex(samples)
L["sex"] = meta.sex.reindex(samples)
L["fib"] = pd.to_numeric(meta.fibrosis_stage.reindex(samples), errors="coerce")
L["nas"] = pd.to_numeric(meta.nas_score.reindex(samples), errors="coerce")
m130 = L.dataset.eq("GSE130970")
L.loc[m130, "group_binary"] = np.where(L.index[m130].map(s130.source_control_status) == "Control", "Control", "Disease")
L.loc[m130, "fib"] = L.index[m130].map(s130.fibrosis_stage).astype(float)
L.loc[m130, "nas"] = L.index[m130].map(s130.nafld_activity_score).astype(float)
L.loc[m130, "sex"] = L.index[m130].map(s130.sex_geo).map({"female": "F", "male": "M"})
g = L.dataset.eq("GSE213621")
L["fib_cat"] = L.fib.map(lambda v: f"F{int(v)}" if pd.notna(v) else "NR")
L.loc[g, "fib_cat"] = L.loc[g, "condition"].map({"Control": "Control", "Fibrosis_F0F1": "F0-F1", "Fibrosis_F2": "F2", "Fibrosis_F3F4": "F3-F4"})
L.loc[g, "fib"] = np.where(L.loc[g, "fib_cat"].eq("F2"), 2.0, np.nan)
L["full_called"] = calls.full_called.reindex(samples)
L["comp"] = L.index.map(cs); L["screen"] = L.index.map(ca)
L["reason"] = np.where(~callable_.reindex(samples).values, "low_sites", "")
bad = {ca[p[0]] for p in failed}
for c, grp in L.groupby("comp"):
    if any(key(a, b) not in stable for a, b in itertools.combinations(grp.index, 2)):
        bad |= set(grp.screen)
L.loc[L.screen.isin(bad), "reason"] = "unresolved_component"
for c, grp in L.groupby("screen"):
    if len(grp) > 1:
        cf = [x for x in ["dataset", "group_binary", "fib_cat", "nas", "sex", "diagnosis_harmonized"] if grp[x].dropna().astype(str).nunique() > 1]
        if cf:
            L.loc[grp.index, "reason"] = "conflict:" + ",".join(cf)
bulk = pd.read_csv(R / "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20260923T231939Z/primary_source_controls/model_design.tsv", sep="\t").sample_id
B = L.loc[bulk].copy()
e = B[B.reason.eq("")].copy(); e["sid"] = e.index
sel = e.sort_values(["comp", "full_called", "sid"], ascending=[True, False, True]).groupby("comp").head(1).index
B["selected"] = B.index.isin(sel)
say(B[B.selected].groupby(["dataset", "group_binary"]).size().unstack().to_string())
say("total", int(B.selected.sum()), "control", int((B.selected & B.group_binary.eq("Control")).sum()), "disease", int((B.selected & B.group_binary.eq("Disease")).sum()))
adopted = set(pd.read_csv(Q / "bulk_primary_sample_ids.tsv", sep="\t").sample_id)
say("identical to adopted 774:", set(sel) == adopted, "mine-only", sorted(set(sel) - adopted), "adopted-only", sorted(adopted - set(sel)))
B["why"] = np.where(B.selected, "selected", np.where(B.reason.ne(""), B.reason, "repeat_not_selected"))
say(pd.crosstab([B.dataset, B.group_binary], B.why).to_string())
gg = B[B.dataset.eq("GSE126848")]
say("GSE126848 full_called by group:\n" + gg.groupby("group_binary").full_called.describe()[["count", "min", "50%", "max"]].to_string())
S = B[B.selected & B.dataset.isin(["GSE130970", "GSE135251", "GSE162694"]) & B.fib.notna()]
say("exact-stage units", len(S), S.fib.astype(int).value_counts().sort_index().to_dict(),
    "F0 control", int((S.fib.eq(0) & S.group_binary.eq("Control")).sum()), "F0 disease", int((S.fib.eq(0) & S.group_binary.eq("Disease")).sum()))
N = B[B.selected & B.dataset.isin(["GSE130970", "GSE135251", "GSE162694"]) & B.nas.notna()]
say("exact-NAS units", len(N), "both stage and NAS", int((N.fib.notna()).sum()))
m213 = B[B.dataset.eq("GSE213621")].groupby("screen").filter(lambda d: len(d) > 1)
r = m213.groupby("screen").agg(n=("fib_cat", "size"), cats=("fib_cat", lambda x: ",".join(sorted(set(map(str, x))))), nsel=("selected", "sum"))
say("GSE213621 multi-library components", len(r), "discordant stage", int(r.cats.str.contains(",").sum()),
    "discordant with any selected", int((r.cats.str.contains(",") & (r.nsel > 0)).sum()),
    "concordant collapsed to exactly one", int((~r.cats.str.contains(",") & r.nsel.eq(1)).sum()), "sizes", r.n.value_counts().to_dict())
