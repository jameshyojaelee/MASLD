import pandas as pd, numpy as np
RUN="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321"
arch=pd.read_csv(f"{RUN}/credible_set_pip_architecture.tsv",sep="\t")
g = pd.read_csv("out_D/D_gated_instances_with_clusters.tsv.gz", sep="\t", dtype={"chr": str})
fixmap = g.set_index("credible_set_uid").pos
arch["chr"]=arch.locus.astype(str).str.split(".").str[0]
arch["rawpos"]=arch.locus.astype(str).str.split(".").str[1].astype(int)
arch["anc5"]=arch.ancestry.replace({"CSA":"SAS"})
arch["gate"]=arch.pip_sum_gate_passed.astype(str).str.lower()=="true"
def ncl(df,col,w):
    n=0
    for c,gg in df.groupby("chr"):
        p=np.sort(gg[col].unique()); n+=1+int(np.sum(np.diff(p)>w))
    return n
target={'ALT':122,'AST':115,'GGT':93,'NAFLD':13,'NASH':4,'PDFF':13}
for gate in [False, True]:
    d = arch[arch.gate] if gate else arch
    d = d.copy(); d["fixpos"]=d.credible_set_uid.map(fixmap)
    for col in ["rawpos"]:
        for w in [2.5e5,5e5,1e6]:
            tt={tr: ncl(x,col,w) for tr,x in d.groupby("trait")}
            ta={(tr,an): ncl(x,col,w) for (tr,an),x in d.groupby(["trait","anc5"])}
            print("gate",gate,col,int(w),tt,sum(tt.values()),"global",ncl(d,col,w))
            if tt==target: print("   MATCH", ta)
