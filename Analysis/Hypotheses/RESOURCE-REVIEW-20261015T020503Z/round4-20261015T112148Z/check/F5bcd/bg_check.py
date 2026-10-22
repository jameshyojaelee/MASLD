# Independent re-derivation of Fig 5D observed/expected peak mass per pair x lineage, from release files.
import numpy as np, pandas as pd, sys
C="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Multimodal_Program_Projection/candidates/"
prim=pd.read_csv(C+"atac-context-v3-coloc-r2-20260924T053244Z/genetics/context/genetic_lineage_context_primary_pairs.tsv",sep="\t")
prim["signal_pair_index"]=prim.signal_pair_index.astype(str)
keys=prim[["gwas_name","ensembl","signal_pair_index"]].drop_duplicates()
print("pairs",len(keys),"rows",len(prim),flush=True)
parts=[]
for ch in pd.read_csv(C+"atac-context-v3-candidate-2026-08-11-r1/genetics/context/variant_liftover_audit.tsv",sep="\t",
    usecols=["gwas_name","ensembl","signal_pair_index","snp_pp_h4","mapping_status","hg38_chrom","hg38_position_1based"],
    dtype=str,chunksize=3_000_000):
    parts.append(ch.merge(keys,on=["gwas_name","ensembl","signal_pair_index"],how="inner"))
v=pd.concat(parts,ignore_index=True); del parts
v["pp"]=v.snp_pp_h4.astype(float)
v["mapped"]=v.mapping_status.eq("mapped")
v["pos"]=pd.to_numeric(v.hg38_position_1based,errors="coerce")
print("variants",len(v),"mapped",v.mapped.sum(),flush=True)
pk=pd.read_csv(C+"atac-context-v3-candidate-2026-08-11-r1/consensus_peak_manifest.tsv",sep="\t",
    usecols=["lineage","chrom","start0","end","supported_gse244832","supported_gse281367"],dtype={"supported_gse244832":str,"supported_gse281367":str})
def flag(sub):
    out=np.zeros(len(v),bool)
    for c,g in sub.groupby("chrom"):
        g=g.sort_values("start0"); s=g.start0.to_numpy(); e=np.maximum.accumulate(g.end.to_numpy())
        m=(v.hg38_chrom.to_numpy()==c)&v.mapped.to_numpy()
        p=v.pos.to_numpy()[m]
        j=np.searchsorted(s,p-1,side="right")-1  # last peak with start0 <= p-1
        ok=j>=0; r=np.zeros(len(p),bool); r[ok]=e[j[ok]]>=p[ok]
        out[np.where(m)[0]]=r
    return out
k=v.gwas_name+"|"+v.ensembl+"|"+v.signal_pair_index
res=[]
for lin in ["hepatocyte","stellate","macrophage","cholangiocyte","t_nk"]:
    L=pk[pk.lineage==lin]
    a=flag(L[L.supported_gse244832.str.upper()=="TRUE"]); b=flag(L[L.supported_gse281367.str.upper()=="TRUE"])
    d=pd.DataFrame({"k":k,"pp":v.pp,"m":v.mapped,"a":a,"b":b})
    d["any"]=d.a|d.b; d["sh"]=d.a&d.b
    for c in ["a","b","any","sh"]: d["pp_"+c]=d.pp*d[c]
    g=d.groupby("k").agg(n_mapped=("m","sum"),n_any=("any","sum"),o244=("pp_a","sum"),o281=("pp_b","sum"),oany=("pp_any","sum"),osh=("pp_sh","sum")).reset_index()
    g["lineage"]=lin; res.append(g)
r=pd.concat(res); r[["gwas_name","ensembl","signal_pair_index"]]=r.k.str.split("|",expand=True)
r["exp_any"]=r.n_any/r.n_mapped
r=prim.merge(r,on=["gwas_name","ensembl","signal_pair_index","lineage"],how="left",validate="1:1")
for o,t in [("oany","any_accessible_mass"),("o244","gse244832_accessible_mass"),("o281","gse281367_accessible_mass"),("osh","shared_accessible_mass")]:
    print("maxabsdiff",o,(r[o]-r[t]).abs().max())
r.drop(columns="k").to_csv(sys.argv[1],sep="\t",index=False)
print("done",len(r))
