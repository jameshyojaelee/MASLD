#!/usr/bin/env python
"""U15 S2I: count of unlabeled genes at best SuSiE PP.H4 >= 0.95 (legend: 144), using the
per-gene best rows from u15_s2a_i.py and the lead-variant class annotation."""
import sys
from pathlib import Path
import pandas as pd
R = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
U = Path(sys.argv[1])
b = pd.read_csv(U / "out_s2ai/S2HI_best_per_gene.tsv", sep="\t")
b = b[b["PP.H4.susie"] > 0.5].copy()
an = pd.read_csv(R / "figures/candidates/fig2-coloc-r2-20260924T053244Z/variant_class/coloc_variant_annotation.csv")
vg = an.drop_duplicates(["variant_key", "ensembl"]).set_index(["variant_key", "ensembl"]).coarse_class
v = an.drop_duplicates("variant_key").set_index("variant_key").coarse_class
cls = [vg.get((t, e), v.get(t, "not annotated")) for t, e in zip(b.top_snp, b.ensembl)]
b["cls"] = cls
print("lead class among 436:", b.cls.value_counts().to_dict())
lab = (b.n_anc09 >= 2) | (b.cls == "coding")
print("labeled", int(lab.sum()), "multi-ancestry", int((b.n_anc09 >= 2).sum()), "coding", int((b.cls == "coding").sum()))
print("unlabeled with best >= 0.95:", int(((b["PP.H4.susie"] >= 0.95) & ~lab).sum()))
print("unlabeled with best > 0.95:", int(((b["PP.H4.susie"] > 0.95) & ~lab).sum()))
b.to_csv(U / "out_s2ai/S2I_labels_rederived.tsv", sep="\t", index=False)
