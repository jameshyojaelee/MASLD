"""audit-H58: count pooled physical clusters that hold EUR or AMR rows (coding and context sets).
Counts only; computes no fraction, difference or ratio."""
import sys, pandas as pd
RUN="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/"
sys.path.insert(0, RUN+"snapshot/scripts/figures")
import fig2_multiplicity_correction as m
a=pd.read_csv(RUN+"multiplicity/corrected_locus_coordinates.tsv",sep="\t",dtype={"locus":str})
p=m.assign_clusters(a.copy())
print("coding: pooled clusters", p.cluster_id.nunique(), "holding EUR/AMR rows", p[p.ancestry.isin(["EUR","AMR"])].cluster_id.nunique())
base=RUN+"snapshot/GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321/"
arch=pd.read_csv(base+"credible_set_pip_architecture.tsv",sep="\t",dtype={"locus":str})
mem=pd.read_csv(base+"credible_set_members_annotated.tsv",sep="\t",low_memory=False,dtype={"locus":str})
arch=arch[arch.trait_scope.isin(m.SCOPES)]; mem=mem[mem.trait_scope.isin(m.SCOPES)]
arch,mem,exc=m.apply_gate(arch,mem)
ps=m.assign_clusters(m.parse_locus(m.build_context_masses(mem),mem))
anc=mem.drop_duplicates("credible_set_uid").set_index("credible_set_uid").ancestry
ps["ancestry"]=ps.credible_set_uid.map(anc)
print("context: pooled clusters", ps.cluster_id.nunique(), "holding EUR/AMR rows", ps[ps.ancestry.isin(["EUR","AMR"])].cluster_id.nunique())
eu=mem[mem.ancestry=="EUR"]
print("EUR member rows", len(eu), "unresolved_liftover rows", int((eu.promoter_context_status=="unresolved_liftover").sum()))
print("EUR studies with unresolved rows:", eu[eu.promoter_context_status=="unresolved_liftover"].groupby(["trait_scope","study"]).size().to_dict())
