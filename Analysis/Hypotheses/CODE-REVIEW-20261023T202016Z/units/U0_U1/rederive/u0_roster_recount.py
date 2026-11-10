"""U0 independent recount of the 774-unit bulk roster.

Written from the Methods text (Units of inference), not from the producer.
Genotype calls are re-read from the merged VCF and KING-robust kinship is
recomputed here with numpy; clusters, label conflicts and the representative
choice are then re-derived and compared with the adopted roster.
"""
import hashlib
import io
import itertools
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

R = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
OUT = Path(sys.argv[1])
OUT.mkdir(parents=True, exist_ok=True)
BCFTOOLS = "/gpfs/commons/home/jameslee/micromamba/envs/allelic_qtl/bin/bcftools"

VCF = R / "RNA-seq/Human/Patient_Cohorts/results/rna_genotypes_restricted/a1-full-20260924T025858Z/vcf/identity_panel.merged.vcf.gz"
SITES = R / "Analysis/MASLD_Model_Benchmark/executions/in-biopsy-allelic-obs-20260926T155651Z/qc/v2-20260929T143503Z/identity_site_annotation.tsv"
A1 = R / "Analysis/MASLD_Model_Benchmark/executions/in-biopsy-allelic-a1-full-20260924T025858Z/identity"
META = R / "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv"
S130 = R / "RNA-seq/Human/Patient_Cohorts/metadata/source/GSE130970/GSE130970_source_diagnosis.tsv"
BULK844 = R / "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20260923T231939Z/primary_source_controls/model_design.tsv"
QUAL = R / "figures/candidates/resource-corrections-20261005-v1/identity_resource_qualified"
ADOPTED_IDS = QUAL / "bulk_primary_sample_ids.tsv"
ADOPTED_ROSTER = QUAL / "bulk_roster.tsv"
ENDPOINT_IDS = QUAL / "endpoint_primary_sample_ids.tsv"

log = open(OUT / "u0_log.txt", "w")


def say(*a):
    msg = " ".join(str(x) for x in a)
    print(msg)
    log.write(msg + "\n")
    log.flush()


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


inputs = [VCF, SITES, A1 / "sample_to_individual.tsv", A1 / "kinship_related_pairs.tsv", META, S130, BULK844,
          ADOPTED_IDS, ADOPTED_ROSTER, ENDPOINT_IDS]
pd.DataFrame({"path": [str(p) for p in inputs], "sha256": [sha(p) for p in inputs]}).to_csv(OUT / "input_sha256.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- sites
ann = pd.read_csv(SITES, sep="\t", keep_default_na=False, na_values=["", "NA"])
ann = ann[ann.error_site.astype(str).str.lower().eq("true") & ann.chrom.isin([f"chr{i}" for i in range(1, 23)])].copy()
ann["half_sha"] = [int(hashlib.sha256(f"{c}:{int(p)}".encode()).hexdigest()[:2], 16) % 2 for c, p in zip(ann.chrom, ann.pos)]
say("autosomal identity sites:", len(ann), "half split sha:", ann.half_sha.value_counts().to_dict(),
    "agreement with annotation 'half' column:", float((ann.half_sha == ann.half).mean()))
regions = OUT / "sites.regions.tsv"
ann[["chrom", "pos"]].to_csv(regions, sep="\t", header=False, index=False)

# ---------------------------------------------------------------- genotypes
samples = subprocess.run([BCFTOOLS, "query", "-l", str(VCF)], capture_output=True, text=True, check=True).stdout.split()
q = subprocess.run([BCFTOOLS, "query", "-R", str(regions), "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT:%DP]\n", str(VCF)],
                   capture_output=True, text=True, check=True).stdout
rows = [line.split("\t") for line in q.rstrip("\n").split("\n")]
say("VCF records returned for sites:", len(rows), "samples:", len(samples))
chrom = np.array([r[0] for r in rows]); pos = np.array([int(r[1]) for r in rows])
alt = np.array([r[3] for r in rows])
keep_rec = np.array([("," not in a) for a in alt])
say("records dropped as multi-allelic:", int((~keep_rec).sum()))
G = np.full((len(rows), len(samples)), -1, dtype=np.int8)
for i, r in enumerate(rows):
    for j, f in enumerate(r[4:]):
        gt, _, dp = f.partition(":")
        if dp in (".", "") or int(dp) < 10:
            continue
        a = gt.replace("|", "/").split("/")
        if len(a) != 2 or "." in a:
            continue  # half calls and missing -> missing
        G[i, j] = int(a[0] != "0") + int(a[1] != "0")
G = G[keep_rec]; chrom = chrom[keep_rec]; pos = pos[keep_rec]
half_map = dict(zip(zip(ann.chrom, ann.pos), ann.half_sha))
half = np.array([half_map.get((c, p), -1) for c, p in zip(chrom, pos)])
say("records with no exact site match (overlap only):", int((half < 0).sum()))
# overlapping records (indel spanning a site) get the half of the site they overlap: keep them in full,
# assign half by nearest site position
if (half < 0).any():
    say("WARNING overlap-only records present; assigned to full only")


def king(Gs):
    called = (Gs >= 0).astype(np.float32)
    het = (Gs == 1).astype(np.float32)
    hom0 = (Gs == 0).astype(np.float32)
    hom2 = (Gs == 2).astype(np.float32)
    nsnp = called.T @ called
    hethet = het.T @ het
    ibs0 = hom0.T @ hom2 + hom2.T @ hom0
    het_i = het.T @ called  # het in i at sites called in j
    denom = het_i + het_i.T
    with np.errstate(divide="ignore", invalid="ignore"):
        kin = (hethet - 2 * ibs0) / denom
    return called.sum(0), nsnp, hethet, ibs0, kin


res = {}
for label, mask in [("full", np.ones(len(G), bool)), ("half0", half == 0), ("half1", half == 1)]:
    res[label] = king(G[mask])
    say(label, "sites:", int(mask.sum()))
idx = {s: i for i, s in enumerate(samples)}
ncalled = pd.DataFrame({"sample_id": samples, "full_called": res["full"][0].astype(int),
                        "half0_called": res["half0"][0].astype(int), "half1_called": res["half1"][0].astype(int)})
ncalled.to_csv(OUT / "my_callable_site_counts.tsv", sep="\t", index=False)

# ---- validate against producer's PLINK tables (comparison only)
for label in ["full", "half0", "half1"]:
    t = pd.read_csv(QUAL / f"{label}_filtered_related_pairs.tsv", sep="\t")
    i = t["#IID1"].map(idx).values; j = t["IID2"].map(idx).values
    _, nsnp, hethet, ibs0, kin = res[label]
    say(f"[{label}] producer pairs={len(t)} max|dKIN|={np.nanmax(np.abs(kin[i, j]-t.KINSHIP.values)):.2e} "
        f"NSNP equal={np.mean(nsnp[i, j]==t.NSNP.values):.4f} IBS0 equal={np.mean(ibs0[i, j]==t.IBS0.values):.4f} "
        f"HETHET equal={np.mean(hethet[i, j]==t.HETHET.values):.4f}")
pc = pd.read_csv(QUAL / "filtered_callable_site_counts.tsv", sep="\t").merge(ncalled, on="sample_id", suffixes=("_prod", "_mine"))
say("callable-site counts identical to producer:", {c: bool((pc[c+"_prod"] == pc[c+"_mine"]).all()) for c in ["full_called", "half0_called", "half1_called"]})

# ---------------------------------------------------------------- pairs and clusters (Methods rule)
n = len(samples)
callable_ = res["full"][0] >= 1000
iu = np.triu_indices(n, 1)
fN, fK = res["full"][1][iu], res["full"][4][iu]
h0N, h0K = res["half0"][1][iu], res["half0"][4][iu]
h1N, h1K = res["half1"][1][iu], res["half1"][4][iu]
both_callable = callable_[iu[0]] & callable_[iu[1]]
hi = (fK >= 0.354)
say("pairs with full KING >= 0.354 (any coverage):", int(hi.sum()),
    "of which cross-cohort:", "computed below")
a1 = pd.read_csv(A1 / "kinship_related_pairs.tsv", sep="\t")
a1_same = {tuple(sorted(p)) for p in zip(a1["#IID1"], a1["IID2"]) if True} if False else \
    {tuple(sorted((x, y))) for x, y, r in zip(a1["#IID1"], a1["IID2"], a1["relation"]) if r == "same_individual"}
cand = {tuple(sorted((samples[a], samples[b]))) for a, b, h, c, nn in zip(iu[0], iu[1], hi, both_callable, fN) if h and c and nn >= 500}
cand_kingonly = set(cand)
cand |= a1_same
say("candidate same-person pairs: KING-panel", len(cand_kingonly), "A1 same_individual", len(a1_same), "union", len(cand))
pairstat = {}
for a, b in cand:
    i, j = idx[a], idx[b]
    ok = (res["full"][1][i, j] >= 500 and res["full"][4][i, j] >= 0.354 and
          res["half0"][1][i, j] >= 250 and res["half0"][4][i, j] >= 0.354 and
          res["half1"][1][i, j] >= 250 and res["half1"][4][i, j] >= 0.354 and callable_[i] and callable_[j])
    pairstat[(a, b)] = ok
stable = {p for p, ok in pairstat.items() if ok}
failed = {p for p, ok in pairstat.items() if not ok}
say("stable pairs:", len(stable), "failed candidate pairs:", len(failed))
prod_checks = pd.read_csv(QUAL / "filtered_identity_pair_checks.tsv", sep="\t")
prod_stable = {tuple(sorted((a, b))) for a, b, s in zip(prod_checks.sample_a, prod_checks.sample_b, prod_checks.stable_filtered_same_identity) if str(s) == "True"}
say("stable pair set equal to producer:", stable == prod_stable, "mine-only", sorted(stable - prod_stable), "prod-only", sorted(prod_stable - stable))


def comps(nodes, edges):
    parent = {s: s for s in nodes}

    def f(s):
        while parent[s] != s:
            parent[s] = parent[parent[s]]
            s = parent[s]
        return s
    for a, b in edges:
        ra, rb = f(a), f(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)
    return {s: f(s) for s in nodes}


comp_stable = comps(samples, stable)
comp_any = comps(samples, stable | failed)

# ---------------------------------------------------------------- labels (independent of the roster)
meta = pd.read_csv(META, keep_default_na=False, na_values=[""])
s130 = pd.read_csv(S130, sep="\t", keep_default_na=False, na_values=["", "NA"])
lab = meta.set_index("sample_id")
lab = lab[~lab.index.duplicated()]
L = pd.DataFrame(index=samples)
L["dataset"] = lab.dataset.reindex(samples)
L["group_binary"] = lab.group_binary.reindex(samples)
L["fib"] = pd.to_numeric(lab.fibrosis_stage.reindex(samples), errors="coerce")
L["nas"] = pd.to_numeric(lab.nas_score.reindex(samples), errors="coerce")
L["sex_reported"] = lab.sex.reindex(samples)
L["diagnosis"] = lab.diagnosis_harmonized.reindex(samples)
L["condition"] = lab.condition.reindex(samples)
s = s130.set_index("run")
in130 = L.dataset.eq("GSE130970")
L.loc[in130, "group_binary"] = np.where(L.index[in130].map(s.source_control_status).astype(str) == "Control", "Control", "Disease")
L.loc[in130, "fib"] = L.index[in130].map(s.fibrosis_stage).astype(float)
L.loc[in130, "nas"] = L.index[in130].map(s.nafld_activity_score).astype(float)
L.loc[in130, "sex_reported"] = L.index[in130].map(s.sex_geo).map({"female": "F", "male": "M"})
g213 = L.dataset.eq("GSE213621")
L["fib_cat"] = L.fib.map(lambda v: f"F{int(v)}" if pd.notna(v) else "NR")
L.loc[g213, "fib_cat"] = L.loc[g213, "condition"].map({"Control": "Control", "Fibrosis_F0F1": "F0-F1", "Fibrosis_F2": "F2", "Fibrosis_F3F4": "F3-F4"})
L.loc[g213, "fib"] = np.where(L.loc[g213, "fib_cat"].eq("F2"), 2.0, np.nan)
say("GSE213621 condition levels:", L.loc[g213, "condition"].value_counts(dropna=False).to_dict())
say("GSE130970 source controls in panel:", int((L.loc[in130, "group_binary"] == "Control").sum()))

# ---------------------------------------------------------------- exclusions
L["full_called"] = res["full"][0].astype(int)
L["comp"] = [comp_stable[x] for x in samples]
L["screen"] = [comp_any[x] for x in samples]
L["reason"] = np.where(~callable_, "low_sites", "")
bad = set()
for p in failed:
    bad.add(comp_any[p[0]])
for c, grp in L.groupby("comp"):
    if any(tuple(sorted((a, b))) not in stable for a, b in itertools.combinations(grp.index, 2)):
        bad.update(grp.screen)
L.loc[L.screen.isin(bad), "reason"] = "unresolved_component"
conflict_cols = ["dataset", "group_binary", "fib_cat", "nas", "sex_reported", "diagnosis"]
conf = {}
for c, grp in L.groupby("screen"):
    if len(grp) < 2:
        continue
    cf = [col for col in conflict_cols if grp[col].dropna().astype(str).nunique() > 1]
    if cf:
        conf[c] = cf
        L.loc[grp.index, "reason"] = "label_conflict:" + ",".join(cf)

# ---------------------------------------------------------------- bulk context (844) and selection
bulk = pd.read_csv(BULK844, sep="\t").sample_id.tolist()
say("bulk context libraries:", len(bulk))
B = L.loc[bulk].copy()
elig = B[B.reason.eq("")].copy()
elig = elig.sort_values(["full_called"], ascending=False, kind="mergesort")
elig["sid"] = elig.index
elig = elig.sort_values(["comp", "full_called", "sid"], ascending=[True, False, True])
sel = elig.groupby("comp").head(1)
B["selected"] = B.index.isin(sel.index)
B["reason2"] = np.where(B.selected, "selected", np.where(B.reason.ne(""), B.reason, "repeat_not_selected"))
B.to_csv(OUT / "my_bulk_roster.tsv", sep="\t")
say("selected units:", int(B.selected.sum()))
tab = B[B.selected].groupby(["dataset", "group_binary"]).size().unstack(fill_value=0)
say(tab.to_string())
say("exclusions by dataset x group x reason:\n" + pd.crosstab([B.dataset, B.group_binary], B.reason2).to_string())

adopted = pd.read_csv(ADOPTED_IDS, sep="\t").sample_id
mine = set(B.index[B.selected])
say("selected set identical to adopted 774:", mine == set(adopted), "mine-only:", sorted(mine - set(adopted)), "adopted-only:", sorted(set(adopted) - mine))

# label agreement with adopted roster
ar = pd.read_csv(ADOPTED_ROSTER, sep="\t", keep_default_na=False, na_values=["", "NA"]).set_index("sample_id")
sel_ids = sorted(mine & set(adopted))
gb_diff = [x for x in sel_ids if ar.loc[x, "group_binary"] != B.loc[x, "group_binary"]]
say("group_binary differences vs adopted roster among selected:", gb_diff)

# ---------------------------------------------------------------- GSE126848 detail (question a)
g = B[B.dataset.eq("GSE126848")]
say("GSE126848 full_called by group (median, min, max):",
    g.groupby("group_binary").full_called.agg(["median", "min", "max", "count"]).to_string())
say("GSE126848 reasons:", pd.crosstab(g.group_binary, g.reason2).to_string())

# ---------------------------------------------------------------- GSE213621 pairs (question c)
g = B[B.dataset.eq("GSE213621")]
multi = g.groupby("screen").filter(lambda d: len(d) > 1)
summ = multi.groupby("screen").apply(lambda d: pd.Series({"n": len(d), "fib_cats": ",".join(sorted(d.fib_cat.astype(str).unique())),
                                                           "groups": ",".join(sorted(d.group_binary.unique())),
                                                           "n_selected": int(d.selected.sum()), "reason": ";".join(sorted(d.reason2.unique()))}))
summ.to_csv(OUT / "gse213621_repeat_components.tsv", sep="\t")
say("GSE213621 repeat components:", len(summ), "discordant-stage:", int(summ.fib_cats.str.contains(",").sum()),
    "discordant with a selected member:", int((summ.fib_cats.str.contains(",") & (summ.n_selected > 0)).sum()),
    "concordant collapsed to one:", int((~summ.fib_cats.str.contains(",") & (summ.n_selected == 1)).sum()))

# ---------------------------------------------------------------- staged units (Fig 3B)
S = B[B.selected & B.dataset.isin(["GSE130970", "GSE135251", "GSE162694"]) & B.fib.notna()]
say("Fig3B exact-stage units:", len(S), S.fib.astype(int).value_counts().sort_index().to_dict())
f0 = S[S.fib.eq(0)]
say("F0 reference: controls", int(f0.group_binary.eq("Control").sum()), "disease", int(f0.group_binary.eq("Disease").sum()))
Snas = B[B.selected & B.dataset.isin(["GSE130970", "GSE135251", "GSE162694"]) & B.nas.notna()]
say("exact-NAS units:", len(Snas), "both exact stage and NAS:", int((B.selected & B.dataset.isin(["GSE130970", "GSE135251", "GSE162694"]) & B.nas.notna() & B.fib.notna()).sum()))

# ---------------------------------------------------------------- Methods identity claims on the 774
selidx = np.array([idx[x] for x in sorted(mine)])
sub = lambda M: M[np.ix_(selidx, selidx)]
fn = sub(res["full"][1]); fk = sub(res["full"][4]); ib = sub(res["full"][3])
tri = np.triu_indices(len(selidx), 1)
say("774: min shared sites over pairs:", int(fn[tri].min()), "max KING:", float(np.nanmax(fk[tri])))
order = np.argsort(-fk[tri])[:3]
ids = np.array(sorted(mine))
for o in order:
    a, b = tri[0][o], tri[1][o]
    say("  top pair", ids[a], ids[b], B.loc[ids[a], "dataset"], B.loc[ids[b], "dataset"], "KING", round(float(fk[a, b]), 4),
        "NSNP", int(fn[a, b]), "IBS0", int(ib[a, b]),
        "halves", round(float(res["half0"][4][idx[ids[a]], idx[ids[b]]]), 3), round(float(res["half1"][4][idx[ids[a]], idx[ids[b]]]), 3))
say("774 pairs with KING >= 0.0884:", int((fk[tri] >= 0.0884).sum()))
ds = L.dataset.values
hi_pairs = [(iu[0][k], iu[1][k]) for k in np.where(fK >= 0.354)[0]]
say("all-panel pairs >= 0.354:", len(hi_pairs), "cross-cohort:", sum(ds[a] != ds[b] for a, b in hi_pairs))
same_ibs0 = [int(res["full"][3][a, b]) for a, b in hi_pairs]
say("max IBS0 among pairs >= 0.354:", max(same_ibs0))

# ---------------------------------------------------------------- endpoint labels (Fig S3B)
E = B[B.selected].copy()
adv = (E.nas >= 5) | (E.fib >= 3) | E.fib_cat.eq("F3-F4")
strict = E.group_binary.eq("Control") & ~(E.nas > 1) & ~(E.fib > 0) & ~E.fib_cat.isin(["F1", "F2", "F3", "F4", "F2", "F3-F4"])
E["endpoint"] = np.where(adv & E.group_binary.eq("Disease"), "Advanced", np.where(strict, "StrictControl", ""))
E = E[E.endpoint.ne("")]
ct = E.groupby(["dataset", "endpoint"]).size().unstack(fill_value=0)
ok = ct[(ct.min(axis=1) >= 5)].index
E = E[E.dataset.isin(ok)]
say("endpoint arms by cohort:\n" + E.groupby(["dataset", "endpoint"]).size().unstack(fill_value=0).to_string())
ep = pd.read_csv(ENDPOINT_IDS, sep="\t").sample_id
say("endpoint set equal to identity endpoint roster:", set(E.index) == set(ep), "mine-only", len(set(E.index) - set(ep)), "roster-only", len(set(ep) - set(E.index)))
E[["dataset", "endpoint"]].to_csv(OUT / "my_endpoint_labels.tsv", sep="\t", index_label="sample_id")
B.loc[B.selected, ["dataset", "group_binary", "sex_reported", "fib", "nas", "fib_cat", "full_called"]].to_csv(OUT / "my_774_labels.tsv", sep="\t", index_label="sample_id")
say("done")
