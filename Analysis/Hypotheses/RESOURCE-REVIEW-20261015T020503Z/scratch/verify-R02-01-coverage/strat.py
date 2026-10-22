import sys, numpy as np, pandas as pd
RUN="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816"
sys.path.insert(0, RUN+"/snapshot/scripts/figures")
import fig2_multiplicity_correction as M
base=RUN+"/snapshot/GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321/"
arch=pd.read_csv(base+"credible_set_pip_architecture.tsv",sep="\t",dtype={"locus":str})
mem=pd.read_csv(base+"credible_set_members_annotated.tsv",sep="\t",low_memory=False,dtype={"locus":str})
arch=arch[arch.trait_scope.isin(M.SCOPES)].copy(); mem=mem[mem.trait_scope.isin(M.SCOPES)].copy()
arch,mem,exc=M.apply_gate(arch,mem)
arch=M.parse_locus(arch,mem)
reg=pd.read_csv(RUN+"/snapshot/GWAS/finemapping/config/gwas_registry.tsv",sep="\t")
arch=arch.merge(reg[["study_name","ld_panel","ancestry"]].rename(columns={"study_name":"study","ancestry":"reg_ancestry"}),on="study",how="left")
print("ancestry col == registry ancestry:", (arch.ancestry==arch.reg_ancestry).all(), "missing ld_panel:", arch.ld_panel.isna().sum())
arch["coding_mass"]=arch[M.CODING_COLS].sum(axis=1); arch["arch_total"]=arch[M.ARCH_COLS].sum(axis=1)
print("\nInstances by scope x ancestry:"); print(pd.crosstab(arch.trait_scope, arch.ancestry, margins=True))
print("\nInstances by scope x ld_panel:"); print(pd.crosstab(arch.trait_scope, arch.ld_panel, margins=True))
print("\nDirect instances by study:"); print(arch[arch.trait_scope==M.SCOPES[0]].groupby(["study","ancestry","ld_panel"]).size())
print("\nEnzyme instances by study:"); print(arch[arch.trait_scope==M.SCOPES[1]].groupby(["study","ancestry","ld_panel"]).size())
def run(frame,label):
    f=M.assign_clusters(frame.copy())
    tab=M.cluster_table(f,"coding_mass","arch_total")
    r=M.analyse(tab,M.fit_rng("coding_mass"))
    keys=["n_clusters_direct","n_clusters_enzyme","n_physical_clusters","n_shared_clusters","fraction_direct","fraction_enzyme","ratio","ratio_boot_lo","ratio_boot_hi","difference","diff_boot_lo","diff_boot_hi","n_boot_finite"]
    print(f"\n[{label}] n_inst direct={ (frame.trait_scope==M.SCOPES[0]).sum()} enzyme={(frame.trait_scope==M.SCOPES[1]).sum()}")
    print("  "+"  ".join(f"{k}={r.get(k)!r}" if not isinstance(r.get(k),float) else f"{k}={r.get(k):.6g}" for k in keys))
    return f,r
full,_=run(arch,"pooled (reproduce headline)")
# direct cluster ancestry composition
d=full[full.trait_scope==M.SCOPES[0]]
print("\nDirect clusters: ancestries present per cluster")
cc=d.groupby("cluster_id").ancestry.apply(lambda s: ",".join(sorted(set(s))))
print(cc.value_counts())
print("Direct clusters containing any EUR instance:", d.groupby("cluster_id").ancestry.apply(lambda s:(s=="EUR").any()).sum())
print("Direct clusters EUR-only:", d.groupby("cluster_id").ancestry.apply(lambda s:(s=="EUR").all()).sum())
run(arch[arch.ancestry=="EUR"],"EUR only")
run(arch[arch.ancestry!="EUR"],"non-EUR only")
for a in sorted(arch.ancestry.unique()):
    run(arch[arch.ancestry==a],f"ancestry={a}")
mvp=arch[arch.study.str.contains("MVP",case=False)]
print("\nMVP studies:", sorted(mvp.study.unique()))
run(mvp[mvp.ancestry!="EUR"],"MVP non-EUR (AFR/AMR/EAS)")
run(mvp,"MVP all ancestries")
for p in sorted(arch.ld_panel.unique()):
    run(arch[arch.ld_panel==p],f"ld_panel={p}")
