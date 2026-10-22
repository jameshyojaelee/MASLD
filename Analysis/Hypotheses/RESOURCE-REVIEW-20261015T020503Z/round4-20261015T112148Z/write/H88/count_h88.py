# H88 planning counts: class sizes among measurable programs. Reads edge metadata and gate measurability only; no spatial statistic.
import pandas as pd, numpy as np, itertools, math
R="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
T4=f"{R}/RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-PROGRAM-MAP-v9/systems/T4_interaction"
G=f"{R}/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/GAP-tissue/out"
dmap=pd.read_csv(f"{R}/Analysis/Hypotheses/V4-20261008T015029Z/hmsma_gsa_donor_map.tsv",sep="\t")
print("donor map groups", dmap.group.value_counts().to_dict(), "| sex x group", pd.crosstab(dmap.group,dmap.sex).to_dict())
prim=set(dmap[dmap.group=="primary"].array_id)
m=pd.read_csv(f"{G}/program_measurability.tsv",sep="\t"); mp=m[m.array_id.isin(prim)]
ok=mp.groupby("program").measurable_080.all(); ok=set(ok[ok].index)
mn8=mp[mp.program.isin(ok)].groupby("program").n_measured.min()
print("measurable all 32:",len(ok),"| min n_measured among them",mn8.min(),"| programs with <8 measured in any donor",int((mn8<8).sum()))
e=pd.read_csv(f"{T4}/program_interaction_edges.tsv.gz",sep="\t")
c=pd.read_csv(f"{T4}/coordinated_edges.tsv",sep="\t")
ck=set(map(frozenset,zip(c.program_a,c.program_b)))
print("edge rows",len(e),"coordinated rows",len(c),"coord in edge table",sum(frozenset(x) in ck for x in zip(e.program_a,e.program_b)))
e["coord"]=[frozenset(x) in ck for x in zip(e.program_a,e.program_b)]
e["meas"]=e.program_a.isin(ok)&e.program_b.isin(ok)
e["lp"]=["-".join(sorted([a,b])) for a,b in zip(e.cell_type_a,e.cell_type_b)]
x=e[e.meas & ~e.same_cell_type]
print("cross-lineage measurable pairs",len(x),"| coord",int(x.coord.sum()),"| entangled",int(x.membership_entangled.sum()),"| n_shared>0",int((x.n_shared>0).sum()))
print(pd.crosstab([x.coord],[x.n_shared>0,x.membership_entangled]))
pop=x[(x.n_shared==0)&(~x.membership_entangled)]
print("POPULATION gene-disjoint, not entangled:",len(pop),"coord",int(pop.coord.sum()),"noncoord",int((~pop.coord).sum()))
print("coord with r_disjoint<0:",int((pop.coord&(pop.r_disjoint<0)).sum()),"| coord r_partial<=0:",int((pop.coord&(pop.r_partial<=0)).sum()))
print("noncoord r_disjoint>0:",int((~pop.coord&(pop.r_disjoint>0)).sum()),"| noncoord r_partial>0:",int((~pop.coord&(pop.r_partial>0)).sum()))
sm=pop[pop.r_disjoint>0]
print("SIGN-MATCHED population (r_disjoint>0):",len(sm),"coord",int(sm.coord.sum()),"noncoord",int((~sm.coord).sum()))
t=pd.crosstab(pop.lp,pop.coord); t2=pd.crosstab(sm.lp,sm.coord); print(pd.concat([t,t2],axis=1,keys=["all","signmatched"]))
print("strata with both classes: all",int(((t>0).all(1)).sum()),"of",len(t),"| sign-matched",int(((t2>0).all(1)).sum()),"of",len(t2))
lin=pd.Series({p:p.split("_")[1] for p in ok}).value_counts(); print("measurable programs by lineage",lin.to_dict())
DUCT="hotspot_hepatocytes_48f39dd4d817a10e"; COLL="hotspot_fibroblasts_98c7913eef5f9315"
print("PC2 pair in population:", any((pop.program_a.isin([DUCT,COLL]))&(pop.program_b.isin([DUCT,COLL]))), "| DUCT measurable", DUCT in ok, "COLL measurable", COLL in ok)
print("edges involving measurable programs per program (coord degree) max", pop[pop.coord].melt(value_vars=["program_a","program_b"]).value.value_counts().max())
