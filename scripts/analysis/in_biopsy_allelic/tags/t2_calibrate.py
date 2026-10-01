#!/usr/bin/env python3
"""Model A pivot, step T2: calibrate the het-calling rule on Geuvadis RNA vs 1kGP DNA.

Input: per-chromosome allele-depth tables at exonic 1kGP SNVs polymorphic among the
30 Geuvadis individuals (columns chrom, pos, ref, alt, then per individual
"refreads,altreads"), and the DNA genotypes of the same individuals.
Rule grid: DP >= d, minor-allele reads >= k, minor-allele fraction >= f.
Choice (fixed before looking at MASLD data): the rule with the highest het
sensitivity among rules whose false-discovery proportion (DNA-hom sites called
het / all sites called het) is <= 0.02; ties go to the larger d, then larger k.
Also reported per rule: sensitivity by RNA allelic-fraction bin and by depth.
"""
import itertools
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

work = Path(sys.argv[1])
recs = []
for f in sorted(work.glob("chr*.ad_vs_dna.tsv.gz")):
    recs.append(pd.read_csv(f, sep="\t"))
x = pd.concat(recs, ignore_index=True)  # one row per site x individual: ref, alt, dna (0/1/2)
x["dp"] = x["ref_reads"] + x["alt_reads"]
x["minor"] = np.minimum(x["ref_reads"], x["alt_reads"])
x["mfrac"] = x["minor"] / x["dp"].where(x["dp"] > 0)
res = []
for d, k, fr in itertools.product((10, 20, 30), (1, 2, 3), (0.01, 0.02, 0.05, 0.1)):
    m = x["dp"] >= d
    call = m & (x["minor"] >= k) & (x["mfrac"] >= fr)
    het = x["dna"] == 1
    called = int(call.sum())
    fdp = float((call & ~het).sum() / max(called, 1))
    sens = float((call & het & m).sum() / max((het & m).sum(), 1))
    res.append({"min_dp": d, "min_minor_reads": k, "min_minor_frac": fr, "het_sensitivity": sens,
                "false_discovery_prop": fdp, "n_dna_het_at_depth": int((het & m).sum()), "n_called": called})
r = pd.DataFrame(res)
ok = r[r["false_discovery_prop"] <= 0.02]
best = ok.sort_values(["het_sensitivity", "min_dp", "min_minor_reads"], ascending=False).iloc[0].to_dict() if len(ok) else None
out = {"grid": r.to_dict("records"), "chosen_rule": best}
if best:
    m = x["dp"] >= best["min_dp"]
    call = m & (x["minor"] >= best["min_minor_reads"]) & (x["mfrac"] >= best["min_minor_frac"])
    h = x[(x["dna"] == 1) & m].assign(called=call[(x["dna"] == 1) & m])
    h["rna_af_bin"] = pd.cut((h["alt_reads"] / h["dp"] - 0.5).abs(), [0, 0.1, 0.2, 0.3, 0.4, 0.45, 0.51], include_lowest=True)
    out["chosen_sensitivity_by_rna_imbalance"] = {str(k): [float(v.called.mean()), int(len(v))] for k, v in h.groupby("rna_af_bin", observed=True)}
    h["dp_bin"] = pd.cut(h["dp"], [0, 20, 30, 50, 100, 10**9])
    out["chosen_sensitivity_by_depth"] = {str(k): [float(v.called.mean()), int(len(v))] for k, v in h.groupby("dp_bin", observed=True)}
    homs = x[(x["dna"] != 1) & m]
    out["chosen_median_ref_fraction_at_true_hets"] = float((h["ref_reads"] / h["dp"]).median())
    out["n_sites_x_individuals"] = int(len(x)); out["n_dna_hom_at_depth"] = int(len(homs))
(work.parent / "t2_calibration.json").write_text(json.dumps(out, indent=2, default=str))
print(json.dumps({k: v for k, v in out.items() if k != "grid"}, indent=2, default=str))
