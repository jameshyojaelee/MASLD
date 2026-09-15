#!/usr/bin/env python
"""A1 supplement: does the element-level activity convention change A1.2?

The spec defines activity per construct. The main script aggregates the two constructs
of an element to the element mean so that activity and allele effect share one
inferential unit (element x context). This script reports the construct-level
alternative, in which the ref and alt constructs are stacked as separate units, so the
A1.2 comparison can be read under either convention.

Usage: python a1_activity_convention_check.py <output_dir>
"""
from __future__ import annotations

import os
import sys
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=RuntimeWarning)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from a1_mpra_reliability import (  # noqa: E402
    CONTEXTS, REPS, build_matrices, half_means, load, reliability_from_halves,
)


def main(outdir):
    tab = os.path.join(outdir, "tables")
    os.makedirs(tab, exist_ok=True)
    el, ro, paired, grp, g2b, o, qc = load()
    elements = sorted(paired.element_id)
    M = build_matrices(o, elements)
    rows = []
    for ctx in CONTEXTS:
        aref = M[(ctx, "ref", "act")].values
        aalt = M[(ctx, "alt", "act")].values
        for name, Mat in (("activity_element_mean_of_alleles", (aref + aalt) / 2.0),
                          ("activity_construct_level_stacked", np.vstack([aref, aalt])),
                          ("activity_ref_construct_only", aref),
                          ("activity_alt_construct_only", aalt),
                          ("allele_effect", aalt - aref)):
            H = np.stack([half_means(Mat, p) for p in range(3)], axis=1)
            ok = np.isfinite(H).all(axis=(1, 2))
            st = reliability_from_halves(H[ok])
            rows.append(dict(context=ctx, quantity=name, n_units=int(ok.sum()),
                             pearson=st["pearson"], spearman=st["spearman"],
                             projected_4rep=st["projected_4rep"]))
    R = pd.DataFrame(rows)
    R.to_csv(f"{tab}/activity_convention_sensitivity.tsv", sep="\t", index=False)
    print(R.round(4).to_string(index=False))
    piv = R.pivot(index="context", columns="quantity", values="pearson")
    verdicts = {c: bool(min(piv.loc[c, q] for q in piv.columns if q.startswith("activity"))
                        > piv.loc[c, "allele_effect"]) for c in CONTEXTS}
    print("\nA1.2 holds under every activity convention, per context:", verdicts)
    with open(f"{tab}/activity_convention_sensitivity.note.txt", "w") as fh:
        fh.write(str(verdicts) + "\n")


if __name__ == "__main__":
    main(sys.argv[1])
