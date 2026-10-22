import pandas as pd, numpy as np
P="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
pairs=pd.read_csv(P+"Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-coloc-r2-20260924T053244Z/genetics/context/genetic_lineage_context_primary_pairs.tsv",sep="\t")
print(len(pairs), pairs.groupby(["gwas_name","gene"]).ngroups, pairs.drop_duplicates(["gwas_name","gene"]).trait_class.value_counts().to_dict())
print(pairs.lineage.value_counts().to_dict())
# relation of state to masses
pairs["m1"]=pairs.gse244832_accessible_mass; pairs["m2"]=pairs.gse281367_accessible_mass
pairs["max1"]=pairs[["m1","m2"]].max(axis=1)
print(pairs.groupby("evidence_state").agg(n=("gene","size"),min_max=("max1","min"),max_max=("max1","max"),min_any=("any_accessible_mass","min"),max_any=("any_accessible_mass","max"),both=("m1",lambda s: ((s>=.5)&(pairs.loc[s.index,"m2"]>=.5)).sum())))
pairs["open_mine"]=pairs.max1>=0.5
print(pd.crosstab(pairs.evidence_state,pairs.open_mine))
reg=pd.read_csv(P+"figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014/program_registry_v2.tsv",sep="\t")
mem=pd.read_csv(P+"Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot/program_membership_v2.tsv",sep="\t")
# membership hash consistent with registry uid
mem["uid2"]="hotspot_"+mem.cell_type+"_"+mem.membership_sha256.str[:16]
print("uid match", (mem.uid2==mem.program_uid).all(), "reg uids in mem", reg.program_uid.isin(mem.program_uid).all(), mem.program_uid.nunique())
print(mem.mapped_symbol_status.value_counts().to_dict())
LIN={"hepatocytes":"hepatocyte","cholangiocytes":"cholangiocyte","macrophages":"macrophage","fibroblasts":"stellate","tcells":"t_nk"}
openp=pairs[pairs.open_mine]
rows=[]
for r in reg.itertuples():
    m=mem[mem.program_uid==r.program_uid]
    genes_all=set(m.mapped_symbol.dropna())
    genes=set(m.loc[m.mapped_symbol_status=="gencode_v49_unique_symbol_confirmed","mapped_symbol"])
    o=openp[(openp.lineage==LIN[r.cell_type])&openp.gene.isin(genes)]
    oa=openp[(openp.lineage==LIN[r.cell_type])&openp.gene.isin(genes_all)]
    rows.append(dict(program_uid=r.program_uid,ct=r.cell_type,module=r.module,name=r.module_name,n_open=o.gene.nunique(),n_open_allsym=oa.gene.nunique(),
      direct=";".join(sorted(o[o.trait_class=="direct_MASLD"].gene.unique())),genes=";".join(sorted(o.gene.unique())),
      enzonly=len(set(o[o.trait_class=="liver_enzyme"].gene)-set(o[o.trait_class=="direct_MASLD"].gene))))
g=pd.DataFrame(rows)
print("programs with open gene:",(g.n_open>0).sum(),"| with all symbols:",(g.n_open_allsym>0).sum())
print("programs with direct open:",(g.direct!="").sum(), g.loc[g.direct!="",["ct","module","name","direct"]].to_string())
print("programs with only-enzyme open genes:",((g.n_open>0)&(g.direct=="")).sum())
g.to_csv(P+"Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/recheck/F5ef/out/my_genetic_117.tsv",sep="\t",index=False)
ws=pd.read_csv(P+"Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/ws/F5ef/out/genetic_context_117.tsv",sep="\t")
mm=g.merge(ws[["program_uid","n_open_genes","open_genes"]],on="program_uid")
print("mismatch with ws:",(mm.n_open!=mm.n_open_genes).sum())
print(mm[mm.n_open!=mm.n_open_allsym][["name","n_open","n_open_allsym"]])
