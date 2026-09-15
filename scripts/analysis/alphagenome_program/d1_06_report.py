#!/usr/bin/env python3
"""D1 final assembly: the matched mapping-bias filter and one summary object.

The prespecified WASP-style site rules are read-level "did any read move" rules. The amendment-01
control arm measures what a single base change does to the same reads when the changed base is NOT
the allele, so the allele arm's excess over that floor is the only part attributable to the allele.
This script adds the matched filter built from that comparison and recomputes the site-level
imbalance under it, through the same site table every other number here comes from.

No model is fitted. Writes tables/mapping_bias_matched.tsv and tables/d1_summary.json.
"""
from __future__ import annotations

import json
import pathlib
import sys

import numpy as np
import pandas as pd

BLOCK_BP = 1_000_000


def main() -> None:
    out = pathlib.Path(sys.argv[1]).resolve()
    t = out / "tables"
    bias = pd.read_csv(t / "mapping_bias_sites.tsv", sep="\t")
    sites = pd.read_csv(t / "lineage_site_summary.tsv", sep="\t")
    feas = pd.read_csv(t / "lineage_feasibility.tsv", sep="\t")
    feasj = json.load((t / "d1_feasibility.json").open())
    biasj = json.load((t / "mapping_bias.json").open())

    # matched filter: the allele arm must not fail more often than the matched non-allelic arm
    bias["excess_fail"] = bias["fail_fraction"] - bias["control_fail_fraction"]
    bias["pass_matched"] = bias["excess_fail"] <= 0
    bias["pass_matched_strict"] = bias["excess_fail"] <= 0.01
    bias.to_csv(t / "mapping_bias_matched.tsv", sep="\t", index=False)

    rows = []
    for cohort in ("GSE281367", "GSE244832"):
        b = bias[bias["cohort"] == cohort].set_index("uid")
        for lineage in ("ALL_READS", "Hepatocyte", "ALL_LABELLED"):
            s = sites[(sites["scope"] == cohort) & (sites["lineage"] == lineage)].copy()
            if s.empty:
                continue
            for label, col in (("unfiltered", None), ("matched", "pass_matched"),
                               ("matched_strict", "pass_matched_strict"),
                               ("wasp_primary", "pass_primary_any_read"),
                               ("wasp_secondary", "pass_secondary_1pct")):
                g = s if col is None else s[s["uid"].map(b[col]).fillna(False).to_numpy(bool)]
                rows.append({"cohort": cohort, "lineage": lineage, "filter": label,
                             "n_sites": int(len(g)), "n_blocks": int(g["block"].nunique()) if len(g) else 0,
                             "mean_log2_alt_over_ref": float(g["mean_log2_alt_over_ref"].mean()) if len(g) else float("nan"),
                             "mean_abs_log2": float(g["mean_log2_alt_over_ref"].abs().mean()) if len(g) else float("nan"),
                             "frac_sites_alt_higher": float((g["mean_log2_alt_over_ref"] > 0).mean()) if len(g) else float("nan")})
    eff = pd.DataFrame(rows)
    eff.to_csv(t / "mapping_bias_effect_matched.tsv", sep="\t", index=False)

    pooled = feas[feas["scope"] == "POOLED_30_DONORS"].sort_values("n_sites_passing_filter", ascending=False)
    summary = {
        "package": "d1-lineage-allelic",
        "prespec_sha256": (out / "d1_prespec.sha256").read_text().strip(),
        "amendments": ["AMENDMENT_01_mapping_bias_control.md"],
        "n_targets": feasj["n_targets"], "n_donors": feasj["n_donors"],
        "n_lineage_streams": feasj["n_lineage_streams"],
        "D1.1": feasj["D1.1"], "D1.2": feasj["D1.2"], "D1.3": feasj["D1.3"],
        "diagnostics_triggered": feasj["diagnostics_triggered"],
        "feasibility_pooled": pooled.to_dict(orient="records"),
        "mapping_bias": biasj,
        "matched_filter": {
            "rule": "a site passes when the allele-swap read failure fraction does not exceed the "
                    "matched non-allelic single-base perturbation failure fraction at the same site",
            "sites_passing": {c: int(bias.loc[bias["cohort"] == c, "pass_matched"].sum())
                              for c in ("GSE281367", "GSE244832")},
            "sites_evaluated": {c: int((bias["cohort"] == c).sum()) for c in ("GSE281367", "GSE244832")},
            "median_excess_fail": {c: float(bias.loc[bias["cohort"] == c, "excess_fail"].median())
                                   for c in ("GSE281367", "GSE244832")},
        },
        "imbalance_under_filters": rows,
    }
    json.dump(summary, (t / "d1_summary.json").open("w"), indent=1, default=float)
    print(eff.to_string(index=False))
    print(json.dumps(summary["matched_filter"], indent=1))


if __name__ == "__main__":
    main()
