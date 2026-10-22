"""H58 plan: count operational clusters per ancestry stratum and class (n only).
Uses the adopted producer's clustering function from the job 22218816 snapshot.
Prints no coding fraction, difference or ratio."""
import sys
import pandas as pd
ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
RUN = ROOT + "figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/"
sys.path.insert(0, RUN + "snapshot/scripts/figures")
import fig2_multiplicity_correction as m

a = pd.read_csv(RUN + "multiplicity/corrected_locus_coordinates.tsv", sep="\t", dtype={"locus": str})
D, E = m.SCOPES
pooled = m.assign_clusters(a.copy())

def counts(sub, label):
    c = m.assign_clusters(sub.copy())
    nd = c.loc[c.trait_scope == D, "cluster_id"].nunique()
    ne = c.loc[c.trait_scope == E, "cluster_id"].nunique()
    nphys = c.cluster_id.nunique()
    nshared = nd + ne - nphys
    # nesting: every within-stratum cluster lies inside exactly one pooled cluster
    j = c[["credible_set_uid", "cluster_id"]].merge(
        pooled[["credible_set_uid", "cluster_id"]].rename(columns={"cluster_id": "pooled_id"}),
        on="credible_set_uid")
    nest = j.groupby("cluster_id").pooled_id.nunique().max()
    print(f"{label:<28} inst_d={int((sub.trait_scope==D).sum()):>4} inst_e={int((sub.trait_scope==E).sum()):>5} "
          f"clus_d={nd:>3} clus_e={ne:>4} physical={nphys:>4} shared={nshared:>3} max_pooled_per_stratum_cluster={nest}")

counts(a, "pooled")
counts(a[a.ancestry == "EUR"], "EUR")
counts(a[a.ancestry != "EUR"], "nonEUR (all pooled)")
for anc in ["AFR", "AMR", "EAS", "SAS"]:
    counts(a[a.ancestry == anc], anc)
counts(a[a.study.str.startswith("MVP_")], "MVP only")
# EUR direct study list (for leave-one-study-out)
print(a[(a.ancestry == "EUR")].groupby(["trait_scope", "study"]).size().to_string())
