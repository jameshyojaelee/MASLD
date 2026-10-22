import sys, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/verify-R02-02-reproduce/snap/scripts/figures')
import numpy as np, pandas as pd
import fig2_multiplicity_correction as F
base='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321/'
arch=pd.read_csv(base+'credible_set_pip_architecture.tsv',sep='\t',dtype={'locus':str})
mem=pd.read_csv(base+'credible_set_members_annotated.tsv',sep='\t',low_memory=False,dtype={'locus':str})
arch=arch[arch.trait_scope.isin(F.SCOPES)].copy(); mem=mem[mem.trait_scope.isin(F.SCOPES)].copy()
arch,mem,exc=F.apply_gate(arch,mem)
def run(members,label,cols=('promoter_mass','abc_mass','accessible_mass','resolved_no_context_mass','unresolved_context_mass'),den='noncoding_pip_mass',extra=None):
    ps=F.build_context_masses(members)
    if extra is not None: ps=extra(ps,members)
    ps=F.assign_clusters(F.parse_locus(ps,members))
    out=[]
    for c in cols:
        r=F.analyse(F.cluster_table(ps,c,den),F.fit_rng(c))
        out.append((label,c,r.get('n_clusters_direct'),r.get('n_clusters_enzyme'),round(r.get('fraction_direct',np.nan),6),round(r.get('fraction_enzyme',np.nan),6),round(r.get('difference',np.nan),6),round(r.get('diff_boot_lo',np.nan),6),round(r.get('diff_boot_hi',np.nan),6)))
    return out
rows=[]
rows+=run(mem,'published_all')
anc=mem.groupby('credible_set_uid').ancestry.first()
rows+=run(mem[mem.ancestry=='EUR'],'EUR_only')
rows+=run(mem[~mem.study.str.startswith('MVP')],'no_MVP')
rows+=run(mem[mem.ancestry!='AMR'],'no_AMR')
# conditional on resolved: denominator = resolved noncoding mass
def add_res(ps,members):
    m=members.copy()
    cons=m.consequence.astype(str)
    is_utr=cons.str.contains(r"(?:^|,)[35]_prime_UTR_variant(?:$|,)",regex=True); is_syn=cons.str.contains(r"(?:^|,)synonymous_variant(?:$|,)",regex=True)
    nc=m[(m.consequence_category=='other_noncoding')|((m.consequence_category=='synonymous_or_utr')&is_utr&~is_syn)].copy()
    res=(nc.promoter_context_status=='resolved')&(nc.lineage_accessibility_status=='resolved')
    nc['res_nc']=nc.normalized_susie_pip.astype(float)*res
    nc['abc_res']=nc.normalized_susie_pip.astype(float)*res*(nc.abc_enhancer_overlap.astype(str).str.lower()=='true')
    s=nc.groupby('credible_set_uid')[['res_nc','abc_res']].sum()
    ps=ps.merge(s,left_on='credible_set_uid',right_index=True,how='left')
    return ps
rows+=run(mem,'resolved_only_denominator',cols=('promoter_mass','abc_res','accessible_mass','resolved_no_context_mass'),den='res_nc',extra=add_res)
df=pd.DataFrame(rows,columns=['subset','context','n_direct_clusters','n_enzyme_clusters','frac_direct','frac_enzyme','difference','diff_lo','diff_hi'])
pd.set_option('display.width',220)
print(df.to_string())
df.to_csv('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/verify-R02-02-reproduce/context_subsets.tsv',sep='\t',index=False)
