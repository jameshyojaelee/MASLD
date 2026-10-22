# H88 audit: metadata-only checks of (a) hub structure of coordinated labels vs program size/detection and
# (b) how many measured genes remain in gene-sharing pairs after shared genes are left out. No spatial statistic.
import pandas as pd, numpy as np, json
from scipy.stats import spearmanr
R="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
T4=f"{R}/RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-PROGRAM-MAP-v9/systems/T4_interaction"
G=f"{R}/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/GAP-tissue/out"
MEM=f"{R}/Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_membership_v2.tsv"
dm=pd.read_csv(f"{R}/Analysis/Hypotheses/V4-20261008T015029Z/hmsma_gsa_donor_map.tsv",sep="\t"); prim=set(dm[dm.group=="primary"].array_id)
m=pd.read_csv(f"{G}/program_measurability.tsv",sep="\t"); mp=m[m.array_id.isin(prim)]
okS=mp.groupby("program").measurable_080.all(); ok=sorted(okS[okS].index)
e=pd.read_csv(f"{T4}/program_interaction_edges.tsv.gz",sep="\t"); c=pd.read_csv(f"{T4}/coordinated_edges.tsv",sep="\t")
key=lambda a,b: tuple(sorted((a,b))); ck={key(a,b) for a,b in zip(c.program_a,c.program_b)}
e["coord"]=[key(a,b) in ck for a,b in zip(e.program_a,e.program_b)]
x=e[e.program_a.isin(ok)&e.program_b.isin(ok)&~e.same_cell_type]
P=x[(x.n_shared==0)&~x.membership_entangled]
deg=pd.concat([P[P.coord].program_a,P[P.coord].program_b]).value_counts().reindex(ok).fillna(0)
tot=pd.concat([P.program_a,P.program_b]).value_counts().reindex(ok)
frac=deg/tot
nmed=mp[mp.program.isin(ok)].groupby("program").n_measured.median().reindex(ok)
ngen=mp[mp.program.isin(ok)].groupby("program").n_genes.first().reindex(ok)
out={}
out["top5_programs_share_of_coord_endpoints_in_P"]=float(deg.sort_values(ascending=False).head(5).sum()/deg.sum())
out["programs_with_coord_fraction_ge_0.5"]=int((frac>=0.5).sum()); out["programs_with_zero_coord_in_P"]=int((deg==0).sum())
rho,p=spearmanr(frac,nmed); out["spearman_coordfrac_vs_median_n_measured"]=[float(rho),float(p)]
rho,p=spearmanr(frac,ngen); out["spearman_coordfrac_vs_n_genes"]=[float(rho),float(p)]
# leave-shared-out remaining measured genes, using actual shared mapped symbols and per-array measured sets is not on disk;
# approximate with membership symbols: shared = intersection of mapped_symbol sets.
mem=pd.read_csv(MEM,sep="\t"); S={p:set(g.mapped_symbol) for p,g in mem[mem.program_uid.isin(ok)].groupby("program_uid")}
sh=x[x.n_shared>0].copy()
sh["n_shared_symbols"]=[len(S[a]&S[b]) for a,b in zip(sh.program_a,sh.program_b)]
nm_min=mp[mp.program.isin(ok)].groupby("program").n_measured.min()
sh["min_left_a_upper"]=[nm_min[a]-len(S[a]&S[b]) for a,b in zip(sh.program_a,sh.program_b)]
sh["min_left_b_upper"]=[nm_min[b]-len(S[a]&S[b]) for a,b in zip(sh.program_a,sh.program_b)]
sh["min_left"]=sh[["min_left_a_upper","min_left_b_upper"]].min(1)
out["gene_sharing_pairs"]=len(sh); out["entangled"]=int(sh.membership_entangled.sum())
out["n_shared_symbols_max"]=int(sh.n_shared_symbols.max()); out["n_shared_col_max"]=int(sh.n_shared.max())
out["pairs_possibly_below_8_after_leaveout (lower bound on genes left = min n_measured minus shared symbols)"]=int((sh.min_left<8).sum())
out["of_which_entangled"]=int(((sh.min_left<8)&sh.membership_entangled).sum())
out["of_which_in_F2c"]=int(((sh.min_left<8)&~sh.membership_entangled).sum())
print(json.dumps(out,indent=1))
