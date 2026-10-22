#!/usr/bin/env python
"""CGB4 inputs: gene-property covariate file and 43 system gene sets (ENSGID).
Universe: genes with continuum and disease statistics, LOEUF and >=1 system
membership; the continuum signature genes are excluded."""
import argparse

import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--cg0", required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()
g = pd.read_csv(f"{a.cg0}/genes.tsv", sep="\t")
w = pd.read_csv(f"{a.cg0}/gene_system_weights.tsv.gz", sep="\t")
sig = g.signature_gene.fillna(False).astype(bool)
u = g[~sig & g.cont_z.notna() & g.disease_t.notna() & g.loeuf.notna() & (g.n_systems > 0)].copy()
cov = pd.DataFrame({"GENE": u.g, "P_cont": u.cont_z.abs(), "P_dis": u.disease_t.abs(),
                    "AveExpr": u.AveExpr, "LOEUF": u.loeuf})
cov.to_csv(f"{a.out}/covar.txt", sep="\t", index=False, na_rep="NA")
w = w[w.gene_id_base.isin(set(u.g))]
with open(f"{a.out}/systems.sets", "w") as f:
    for s, d in w.groupby("community_id"):
        f.write(s + " " + " ".join(sorted(d.gene_id_base)) + "\n")
print(f"universe {len(cov)} genes; signature excluded {int(sig.sum())}; "
      f"systems {w.community_id.nunique()}; corr(P_cont,P_dis)={np.corrcoef(cov.P_cont, cov.P_dis)[0,1]:.3f}")
