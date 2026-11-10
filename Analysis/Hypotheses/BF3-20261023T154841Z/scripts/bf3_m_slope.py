#!/usr/bin/env python
"""BF3 post hoc (DEVIATIONS.md 2026-10-23T16:06Z): diagnosis slope by cell type for M, U and F from donor_terms.tsv. Point estimates."""
import sys

import numpy as np
import pandas as pd

PRIMARY, OUT = sys.argv[1], sys.argv[2]
DIAG_ORD = {"Healthy": 0, "Steatosis": 1, "Steatohepatitis": 2, "Cirrhosis": 3}
t = pd.read_csv(f"{PRIMARY}/donor_terms.tsv", sep="\t")
cols = ["S", "COMP", "WITHIN", "INTER"] + [c for c in t.columns if c.startswith(("COMP_", "WITHIN_"))]
rows = []
for s, d in t.groupby("set"):
    X = np.column_stack([d.diagnosis.map(DIAG_ORD), (d.sex == "Male").astype(float),
                         (d.dataset == "GSE202379").astype(float), (d.dataset == "GSE244832").astype(float)])
    rows.append(pd.Series(np.linalg.lstsq(X, d[cols].values, rcond=None)[0][0], index=cols, name=s))
pd.DataFrame(rows).T.round(4).to_csv(f"{OUT}/diagnosis_slope_by_type.tsv", sep="\t")
print(pd.DataFrame(rows).T.round(3))
