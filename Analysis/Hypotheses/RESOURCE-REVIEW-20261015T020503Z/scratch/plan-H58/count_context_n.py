"""H58 plan: n for the noncoding-context family and the liftover status, by stratum.
Counts rows, instances and clusters only; prints no context fraction or difference."""
import sys
import pandas as pd
ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
RUN = ROOT + "figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/"
sys.path.insert(0, RUN + "snapshot/scripts/figures")
import fig2_multiplicity_correction as m
base = RUN + "snapshot/GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321/"
arch = pd.read_csv(base + "credible_set_pip_architecture.tsv", sep="\t", dtype={"locus": str})
mem = pd.read_csv(base + "credible_set_members_annotated.tsv", sep="\t", low_memory=False, dtype={"locus": str})
print("member columns:", list(mem.columns))
arch = arch[arch.trait_scope.isin(m.SCOPES)].copy(); mem = mem[mem.trait_scope.isin(m.SCOPES)].copy()
arch, mem, exc = m.apply_gate(arch, mem)
print("gated member rows", len(mem), "instances", mem.credible_set_uid.nunique(), "unique sites",
      mem[["chromosome", "position"]].drop_duplicates().shape[0])
per_set = m.parse_locus(m.build_context_masses(mem), mem)
anc = mem.drop_duplicates("credible_set_uid").set_index("credible_set_uid").ancestry
per_set["ancestry"] = per_set.credible_set_uid.map(anc)
D, E = m.SCOPES
def counts(sub, label):
    c = m.assign_clusters(sub.copy())
    nd = c.loc[c.trait_scope == D, "cluster_id"].nunique(); ne = c.loc[c.trait_scope == E, "cluster_id"].nunique()
    print(f"context {label:<8} inst_d={int((sub.trait_scope==D).sum()):>4} inst_e={int((sub.trait_scope==E).sum()):>5} clus_d={nd:>3} clus_e={ne:>4} physical={c.cluster_id.nunique()}")
counts(per_set, "pooled")
for a_ in ["EUR", "AMR", "AFR", "EAS", "SAS"]:
    counts(per_set[per_set.ancestry == a_], a_)
# row counts by status (no PIP)
print(pd.crosstab(mem.ancestry, mem.promoter_context_status))
print(pd.crosstab(mem.ancestry, mem.lineage_accessibility_status))
# positive-control sites: row counts only
for site in [("22", 44324727), ("19", 19379549)]:
    s = mem[(mem.chromosome.astype(str) == site[0]) & (mem.position == site[1])]
    print("site", site, "rows", len(s)); print(s.groupby(["trait_scope", "ancestry"]).size().to_string())
