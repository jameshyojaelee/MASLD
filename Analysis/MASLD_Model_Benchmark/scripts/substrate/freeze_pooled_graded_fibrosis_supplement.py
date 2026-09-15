#!/usr/bin/env python
"""Supplement to the pooled graded-fibrosis substrate freeze.

Adds, without touching any frozen file already written:
  1. Detection-floor fractions for the two STAR-variant proportion cohorts, so
     the variant mismatch is measured rather than assumed.
  2. The compositional-closure diagnostic the brief asks for: per-lineage
     Spearman against fibrosis stage computed on the raw share, on the CLR of
     the full 16-part closure, and on the CLR of the eligible subcomposition.
     This is a SUBSTRATE DIAGNOSTIC, not an endpoint.
Then regenerates SHA256SUMS over the whole frozen directory.
"""
import hashlib, json, os, sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

OUT = sys.argv[1]
DECON = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/results"
FLOOR = 1e-4

def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()

rep = json.load(open(os.path.join(OUT, "substrate_freeze_report.json")))
s = pd.read_csv(os.path.join(OUT, "per_sample_substrate.tsv"), sep="\t")
axis = rep["lineages"]["lineage_axis"]
elig = rep["lineages"]["eligible_lineages"]

sup = {"note": "supplement to substrate_freeze_report.json; adds nothing to the frozen membership"}

# integrity of the deposit
sup["deposit_integrity"] = dict(
    n_rows=int(len(s)),
    n_with_all_16_proportions=int(s["proportions_available"].sum()),
    n_missing_proportions=int((~s["proportions_available"]).sum()),
    expression_row_index_is_identity=bool((s["expression_row_index"].to_numpy() == np.arange(len(s))).all()),
    per_cohort_variant={k: v for k, v in s.groupby("cohort")["proportions_variant"].first().items()},
)

# 1. STAR-variant floor fractions
star_cohorts = [g for g, m in rep["lineages"]["proportions_availability"].items()
                if m["deconvolution_input_quantification"] != "kallisto"]
star_frac = {}
for g in star_cohorts:
    pr = pd.read_csv(f"{DECON}/{g}/{g}_bayesprism_proportions.tsv", sep="\t", index_col=0)
    ids = [i for i in s.loc[s.cohort == g, "sample_id"] if i in pr.index]
    if not ids:
        raise RuntimeError(f"{g}: zero-row join")
    star_frac[g] = (pr.loc[ids] > FLOOR).mean(axis=0).round(4).to_dict()
sup["star_variant_fraction_above_floor"] = star_frac
sup["star_variant_would_be_eligible"] = {
    g: sorted([l for l, f in d.items() if f >= 0.50]) for g, d in star_frac.items()}
sup["variant_agreement_with_kallisto_eligible_set"] = {
    g: dict(kallisto_eligible=elig,
            star_eligible=sorted([l for l, f in d.items() if f >= 0.50]),
            overlap=sorted(set(elig) & {l for l, f in d.items() if f >= 0.50}))
    for g, d in star_frac.items()}

# 2. Closure diagnostic
diag = {}
for g, sub in s.groupby("cohort"):
    y = sub["fibrosis_grade"].to_numpy(dtype=float)
    d = {}
    for lin in elig:
        raw = sub[f"prop__{lin}"].to_numpy(dtype=float)
        c16 = sub[f"clr16__{lin}"].to_numpy(dtype=float)
        csb = sub[f"clrsub__{lin}"].to_numpy(dtype=float)
        d[lin] = dict(
            rho_raw_share=round(float(spearmanr(raw, y).statistic), 4),
            rho_clr_full16=round(float(spearmanr(c16, y).statistic), 4),
            rho_clr_eligible_subcomposition=round(float(spearmanr(csb, y).statistic), 4),
        )
    diag[g] = d
sup["closure_diagnostic"] = dict(
    status="substrate_diagnostic_not_an_endpoint",
    why="Cell-type shares are a closed composition. A share can move because its lineage moved or because another did. Reporting the raw share, the full-16 CLR and the eligible-subcomposition CLR side by side makes a closure-induced sign flip visible instead of silent.",
    per_cohort=diag,
)
flips = []
for g, d in diag.items():
    for lin, v in d.items():
        if np.sign(v["rho_clr_full16"]) != np.sign(v["rho_clr_eligible_subcomposition"]):
            flips.append(dict(cohort=g, lineage=lin, **v))
sup["closure_sign_flips_full16_vs_subcomposition"] = flips

with open(os.path.join(OUT, "substrate_freeze_supplement.json"), "w") as fh:
    json.dump(sup, fh, indent=2)

dig = {}
for fn in sorted(os.listdir(OUT)):
    p = os.path.join(OUT, fn)
    if os.path.isfile(p) and fn != "SHA256SUMS":
        dig[fn] = sha256(p)
with open(os.path.join(OUT, "SHA256SUMS"), "w") as fh:
    for fn, d in dig.items():
        fh.write(f"{d}  {fn}\n")
print(json.dumps(sup, indent=2))
print("\n=== SHA256SUMS ===")
print(open(os.path.join(OUT, "SHA256SUMS")).read())
