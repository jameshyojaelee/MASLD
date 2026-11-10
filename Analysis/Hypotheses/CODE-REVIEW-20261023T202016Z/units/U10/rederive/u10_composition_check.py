"""Login-node table check (<1 s): composition-adjusted call retention for the 18 strict maxT systems
(S4U legend: '14 lose the call in at least one composition-adjusted fit')."""
import pandas as pd
R = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
t = pd.read_csv(f"{R}/figures/candidates/continuum-rare-type-arms-20261012T020132Z/arms/meta_effects_all_arms.tsv", sep="\t")
s = t[t.family == "systems"]
mem = pd.read_csv(f"{R}/RNA-seq/results/histology_anchored_continuum/molecular_systems_inference/"
                  "hac-system-inference-corrected774-20261009T171017Z/tables/system_inference_membership.tsv", sep="\t")
m18 = sorted(mem.loc[mem.permutation_maxT_supported_all_four, "community_id"])
print("strict maxT systems", len(m18))
print("units per arm:", s.groupby(["arm", "fraction_set"]).n.first().to_dict())
comp = s[(s.arm != "U_no_composition") & (s.n == 367)]
k = comp[comp.feature_id.isin(m18)].groupby("feature_id").agg(n_fits=("call", "size"), n_called=("call", "sum"))
print(k.to_string())
print("called in all fits", int((k.n_called == k.n_fits).sum()), "| lose in >=1", int((k.n_called < k.n_fits).sum()),
      "| called in <=5", int((k.n_called <= 5).sum()), "| none", int((k.n_called == 0).sum()))
print("unadjusted (U_no_composition) calls among 43 systems:", int(s[s.arm == "U_no_composition"].call.sum()))
