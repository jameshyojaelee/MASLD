"""Critic check: do the 521 Figure 7B/7D graded-fibrosis/routing samples contain repeat libraries
from one genetic identity, and does any identity span two cohorts (which would leak across the
leave-one-cohort-out folds)?"""
import pandas as pd
EX = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/"
L = pd.read_csv(EX + "pooled_graded_fibrosis_ladder-20260828T114911Z/predictions_per_sample.tsv", sep="\t")
I = pd.read_csv(EX + "in-biopsy-allelic-a1-full-20260924T025858Z/identity/sample_to_individual.tsv", sep="\t")
K = pd.read_csv(EX + "in-biopsy-allelic-a1-full-20260924T025858Z/identity/kinship_related_pairs.tsv", sep="\t")
print("kinship cols:", list(K.columns))
m = L.merge(I[["run", "individual_id", "passes_min_sites"]], left_on="sample_id", right_on="run", how="left")
print("521 samples; in identity table:", m.run.notna().sum(), "; with individual_id:", m.individual_id.notna().sum())
print(m.groupby("cohort").agg(n=("sample_id", "size"), ided=("individual_id", lambda s: s.notna().sum()),
                              uniq=("individual_id", lambda s: s.dropna().nunique())))
g = m.dropna(subset=["individual_id"]).groupby("individual_id")
multi = g.filter(lambda d: len(d) > 1)
print("identities with >1 library among the 521:", multi.individual_id.nunique(), "libraries:", len(multi))
cross = g.cohort.nunique(); print("identities spanning >1 cohort:", int((cross > 1).sum()))
if len(multi):
    print(multi.sort_values("individual_id")[["individual_id", "sample_id", "cohort", "fibrosis_grade"]].to_string())
