#!/usr/bin/env python3
"""F5ef: corrected Fig 5F source table (R06-11, R12-03, R12-07, G3-03, G3-07).

Input: the adopted Fig 5F source table (read only) and review job 1731, which re-derived the
matched-gene residual Moran's I for Stromal ECM and Ductular injury under both adjustment designs
(lineage_only = program lineage abundance + 2 QC covariates; all16 = 16 cell2location q05 abundances
+ 2 QC covariates), 9,999 matched-gene draws, per section, donor (GSE192741) and array (Vu).

Changes relative to the adopted table:
  * evidence_state of Vu Ductular injury: source_dependent -> indeterminate. Source dependence moves
    to its own column (donor_join), so it no longer shares a state with supported Stromal ECM.
  * adds k (null draws >= observed), the empirical P = (k+1)/(n+1) and whether P sits at the 1e-4 floor.
  * adds unit-level counts under the plotted all-16 design: units with positive residual I, units with
    positive excess over their own matched-gene null, and the exact one-sided sign-flip P of the mean
    unit excess (2^n patterns). The gene-set null resamples genes with the donors held fixed; the sign
    flip is the donor-level question and its floor at 4 donors is 1/16 = 0.0625. Vu units are arrays,
    not donors, so the Vu sign flip is reported as array-level and not as donor replication.
  * adds the overlay-section value (JBO019, donor H37), which is not the 4-donor collapsed value.
No new statistic is fitted; every number is read from job 1731 and re-checked against the adopted table.
"""
from __future__ import annotations

import hashlib
import itertools
from pathlib import Path

import numpy as np
import pandas as pd

P = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
SRC = P / "figures/main/fig5_molecular_context/panels/data/fig5f_spatial_program_calibration.tsv"
J1731 = P / "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/1731/out"
CAND = P / "figures/candidates/review-r3-fixes-20261015T112148Z/F5ef"
WS = P / "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/ws/F5ef"

src = pd.read_csv(SRC, sep="\t")
prog = pd.read_csv(J1731 / "program_level.tsv", sep="\t")
unit = pd.read_csv(J1731 / "unit_level.tsv", sep="\t")
assert len(src) == 4 and len(prog) == 8

prog["display_name"] = prog["program_name"].str.replace(r" \(.*\)$", "", regex=True)
unit["display_name"] = unit["program_name"].str.replace(r" \(.*\)$", "", regex=True)
a16 = prog[prog["design"] == "all16"].set_index(["display_name", "dataset"])
lin = prog[prog["design"] == "lineage_only"].set_index(["display_name", "dataset"])


def exact_signflip_one_sided(x: np.ndarray) -> float:
    obs = x.mean()
    means = [np.mean(x * np.array(s)) for s in itertools.product([1, -1], repeat=len(x))]
    return float(np.mean(np.array(means) >= obs - 1e-15))


rows = []
for r in src.itertuples(index=False):
    key = (r.display_name, r.source_dataset_id)
    a = a16.loc[key]
    l = lin.loc[key]
    # the adopted table must equal the all-16 re-derivation
    assert abs(a.residual_moran_i - r.residual_moran_i) < 1e-9, key
    assert abs(a.excess - r.moran_excess) < 1e-6, key
    assert a.matched_set_sha256.startswith(r.matched_set_sha256), key
    u = unit[(unit["display_name"] == r.display_name) & (unit["dataset"] == r.source_dataset_id)
             & (unit["design"] == "all16") & (unit["unit_type"] == "donor_or_array")]
    is_gse = r.source_dataset_id == "GSE192741"
    exc = u["excess"].to_numpy()
    d = r._asdict()
    d["evidence_state_adopted"] = r.evidence_state
    d["within_source_state"] = "supported" if r.supported_within_source else "indeterminate"
    if is_gse:
        d["donor_join"] = "joined: 5 sections collapsed to 4 donors"
        d["evidence_state"] = d["within_source_state"]
    else:
        d["donor_join"] = "not joined: 10 arrays, no public donor key"
        d["evidence_state"] = "source_dependent" if r.supported_within_source else "indeterminate"
    d["adjustment_model"] = ("residualized for 16 cell2location q05 abundances, log library size and "
                             "detected genes")
    d["k_null_ge_obs"] = int(a.k_null_ge_obs)
    d["p_empirical"] = float(a.p)
    d["p_at_floor_1e-4"] = bool(a.k_null_ge_obs == 0)
    d["null_resamples"] = "genes (matched on expression, detection, biotype, lineage correlation); units held fixed"
    d["unit_type"] = "donor" if is_gse else "array"
    d["n_units"] = len(u)
    d["n_units_residual_I_positive"] = int((u["residual_moran_i"] > 0).sum())
    d["n_units_excess_positive"] = int((exc > 0).sum())
    d["unit_excess_values"] = ";".join(f"{i}:{v:.4f}" for i, v in zip(u["unit_id"], exc))
    d["unit_signflip_exact_p_one_sided"] = exact_signflip_one_sided(exc)
    d["unit_signflip_floor"] = 1 / 2 ** len(exc)
    d["lineage_only_residual_moran_i"] = float(l.residual_moran_i)
    d["lineage_only_excess"] = float(l.excess)
    d["lineage_only_k_null_ge_obs"] = int(l.k_null_ge_obs)
    if is_gse and r.display_name == "Stromal ECM":
        s = unit[(unit["display_name"] == "Stromal ECM") & (unit["dataset"] == "GSE192741")
                 & (unit["design"] == "all16") & (unit["unit_type"] == "section") & (unit["unit_id"] == "JBO019")]
        assert len(s) == 1
        d["overlay_section"] = "JBO019 (donor H37)"
        d["overlay_section_residual_moran_i"] = float(s["residual_moran_i"].iloc[0])
        d["overlay_section_excess"] = float(s["excess"].iloc[0])
        d["overlay_section_p"] = float(s["p"].iloc[0])
    else:
        d["overlay_section"] = ""
        d["overlay_section_residual_moran_i"] = np.nan
        d["overlay_section_excess"] = np.nan
        d["overlay_section_p"] = np.nan
    rows.append(d)

out = pd.DataFrame(rows)
assert (out.loc[(out.display_name == "Ductular injury") & (out.source_dataset_id == "Vu_et_al_2025"),
                "evidence_state"] == "indeterminate").all()

# Paired SECM minus DI excess per donor (G3-04): between-program contrast at the donor level
pair = []
for ds in ["GSE192741", "Vu_et_al_2025"]:
    uu = unit[(unit["dataset"] == ds) & (unit["design"] == "all16") & (unit["unit_type"] == "donor_or_array")]
    w = uu.pivot(index="unit_id", columns="display_name", values="excess")
    diff = (w["Stromal ECM"] - w["Ductular injury"]).to_numpy()
    pair.append({"source_dataset_id": ds, "unit_type": "donor" if ds == "GSE192741" else "array",
                 "n_units": len(diff), "n_units_secm_minus_di_positive": int((diff > 0).sum()),
                 "unit_diff_values": ";".join(f"{i}:{v:.4f}" for i, v in zip(w.index, diff)),
                 "signflip_exact_p_one_sided": exact_signflip_one_sided(diff),
                 "signflip_floor": 1 / 2 ** len(diff),
                 "note": ("donor-level between-program contrast" if ds == "GSE192741"
                          else "array-level; arrays are not donors")})
pair = pd.DataFrame(pair)

CAND.mkdir(parents=True, exist_ok=True)
out.to_csv(CAND / "fig5f_spatial_program_calibration_corrected.tsv", sep="\t", index=False)
pair.to_csv(CAND / "fig5f_between_program_unit_contrast.tsv", sep="\t", index=False)
with open(WS / "out" / "fig5f_inputs.sha256", "w") as fh:
    for f in (SRC, J1731 / "program_level.tsv", J1731 / "unit_level.tsv"):
        fh.write(f"{hashlib.sha256(f.read_bytes()).hexdigest()}  {f}\n")
pd.set_option("display.width", 250)
print(out[["display_name", "source_dataset_id", "evidence_state_adopted", "evidence_state", "donor_join",
           "moran_excess", "k_null_ge_obs", "p_empirical", "bh_q", "n_units", "n_units_residual_I_positive",
           "n_units_excess_positive", "unit_signflip_exact_p_one_sided", "overlay_section_excess"]].to_string())
print(pair.to_string())
