import sys, numpy as np, pandas as pd
exec(open(sys.argv[1]).read().split("def run(")[0].split('print("ancestry col')[0])
arch["coding_mass"]=arch[M.CODING_COLS].sum(axis=1); arch["arch_total"]=arch[M.ARCH_COLS].sum(axis=1)
arch["stratum"]=np.where(arch.ancestry=="EUR","EUR","nonEUR")
res={}
for s in ["EUR","nonEUR"]:
    f=M.assign_clusters(arch[arch.stratum==s].copy()); tab=M.cluster_table(f,"coding_mass","arch_total")
    for sc in M.SCOPES:
        t=tab[tab.trait_scope==sc]; res[(s,sc)]=(t.num.sum(),t.den.sum(),len(t))
for k,v in res.items(): print(k, "num=%.3f den=%.3f nclust=%d frac=%.4f"%(v[0],v[1],v[2],v[0]/v[1]))
D=M.SCOPES[0]; E=M.SCOPES[1]
wd={s:res[(s,D)][1] for s in ["EUR","nonEUR"]}; tot=sum(wd.values())
fd=sum(res[(s,D)][0] for s in wd)/tot
fe_std=sum(wd[s]/tot*res[(s,E)][0]/res[(s,E)][1] for s in wd)
print("direct-arm stratum weights", {s:round(w/tot,3) for s,w in wd.items()})
print("direct fraction (strata-summed) %.4f; enzyme fraction standardized to direct ancestry mix %.4f; ratio %.3f"%(fd,fe_std,fd/fe_std))
