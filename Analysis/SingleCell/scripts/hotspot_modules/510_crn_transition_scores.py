#!/usr/bin/env python
# KEY OUTPUT: per-module log-ratio of donor module score across the 4 CRN
# transitions (F0->F1, F1->F2, F2->F3, F3->F4). Lets supplementary figures
# overlay unsupervised Hotspot modules on the multi-step cascade narrative
# (paper_outline.md L17 — replaces single F2-switch framing).
import os
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"))
RES = ROOT / "Analysis/SingleCell/results_gpu_v2/hotspot_modules"
META = ROOT / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_metadata_extended.tsv"
BOOTSTRAP_FLAG = ROOT / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/BOOTSTRAP_FALLBACK_REQUIRED.flag"

scores = pd.read_csv(RES / "donor_scores_all.tsv", sep="\t")
META_COLS = ["sample", "F_stage_augmented", "F_stage_augmented_clean",
             "F_stage_inferred", "F_stage_documented",
             "disease_stage_coarse", "dataset", "exclude_stage_analysis"]
meta_full = pd.read_csv(META, sep="\t")
# Backwards-compat: optional columns added by Phase 1 remediation
for col in ("F_stage_augmented_clean", "exclude_stage_analysis"):
    if col not in meta_full.columns:
        meta_full[col] = meta_full["F_stage_augmented"] if col == "F_stage_augmented_clean" else False
meta = meta_full[META_COLS].copy()

# Protocol remediation: drop excluded donors before F-stage priority chain
meta = meta[~meta["exclude_stage_analysis"].astype(bool)]

# F-stage priority — respects bootstrap gate
# When BOOTSTRAP_FALLBACK_REQUIRED.flag exists (set by 343q bootstrap stability gate),
# F_stage_augmented is unreliable; fall back to documented > inferred (skip augmented).
# Otherwise: documented > F_stage_augmented_clean (NA for excluded) > inferred.
if BOOTSTRAP_FLAG.exists():
    print(f"[510] Bootstrap fallback flag detected at {BOOTSTRAP_FLAG} — skipping F_stage_augmented")
    meta["F_stage"] = meta["F_stage_documented"].combine_first(meta["F_stage_inferred"])
else:
    meta["F_stage"] = meta["F_stage_documented"].combine_first(
        meta["F_stage_augmented_clean"]).combine_first(meta["F_stage_inferred"])
meta = meta.dropna(subset=["F_stage"])
meta["F_stage"] = meta["F_stage"].astype(int)

print(f"[510] After protocol+F-stage filters: {len(meta)} donors")
print(meta["F_stage"].value_counts().sort_index().to_string())

df = scores.merge(meta[["sample", "F_stage", "dataset"]], on="sample", how="inner")

TRANSITIONS = [(0, 1), (1, 2), (2, 3), (3, 4)]
records = []
for (ct, mod), g in df.groupby(["cell_type", "module"]):
    # Module column in donor_scores_all.tsv is namespaced as "<ct>__<int>"
    mod_int = int(str(mod).split("__")[-1])
    for f_from, f_to in TRANSITIONS:
        s_from = g.loc[g.F_stage == f_from, "score"]
        s_to = g.loc[g.F_stage == f_to, "score"]
        if len(s_from) < 3 or len(s_to) < 3:
            continue
        m_from, m_to = s_from.mean(), s_to.mean()
        # Welch t-test approximation via mean difference + pooled std
        pooled = np.sqrt(s_from.var(ddof=1) / len(s_from) +
                         s_to.var(ddof=1) / len(s_to))
        delta = m_to - m_from
        z = delta / pooled if pooled > 0 else 0.0
        records.append({
            "cell_type": ct,
            "module": mod_int,
            "transition": f"F{f_from}_to_F{f_to}",
            "n_from": len(s_from),
            "n_to": len(s_to),
            "mean_from": m_from,
            "mean_to": m_to,
            "delta": delta,
            "z": z,
        })

out = pd.DataFrame(records)
# BH-FDR per transition (since modules are tested in parallel within each)
from scipy.stats import norm
out["p_two_sided"] = 2 * (1 - norm.cdf(np.abs(out["z"])))
out["q"] = np.nan
for trans, sub in out.groupby("transition"):
    p = sub["p_two_sided"].values
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q_final = np.empty_like(q)
    q_final[order] = q
    out.loc[sub.index, "q"] = q_final
out["q"] = out["q"].clip(0, 1)

out_path = RES / "crn_transition_scores.tsv"
out.to_csv(out_path, sep="\t", index=False)
print(f"Wrote {out_path} with {len(out)} (module, transition) rows")
print(out.groupby("transition").size())
