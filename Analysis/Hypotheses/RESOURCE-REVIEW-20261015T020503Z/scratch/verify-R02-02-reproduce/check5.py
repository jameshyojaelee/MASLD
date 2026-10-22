import sys, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/verify-R02-02-reproduce/snap/scripts/figures')
import numpy as np, pandas as pd
import fig2_multiplicity_correction as F
base='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321/'
arch=pd.read_csv(base+'credible_set_pip_architecture.tsv',sep='\t',dtype={'locus':str})
mem=pd.read_csv(base+'credible_set_members_annotated.tsv',sep='\t',low_memory=False,dtype={'locus':str})
arch=arch[arch.trait_scope.isin(F.SCOPES)].copy(); mem=mem[mem.trait_scope.isin(F.SCOPES)].copy()
arch,mem,exc=F.apply_gate(arch,mem)
cons=mem.consequence.astype(str)
is_utr=cons.str.contains(r"(?:^|,)[35]_prime_UTR_variant(?:$|,)",regex=True); is_syn=cons.str.contains(r"(?:^|,)synonymous_variant(?:$|,)",regex=True)
nc=mem[(mem.consequence_category=='other_noncoding')|((mem.consequence_category=='synonymous_or_utr')&is_utr&~is_syn)]
u=nc.promoter_context_status=='unresolved_liftover'
print('noncoding rows (script definition) unresolved row fraction by ancestry (%):')
print((u.groupby(nc.ancestry).mean()*100).round(2).to_dict(), 'n', nc.ancestry.value_counts().to_dict())
anc=mem.groupby('credible_set_uid').ancestry.first(); study=mem.groupby('credible_set_uid').study.first()
ps=F.build_context_masses(mem)
ps['ancestry']=ps.credible_set_uid.map(anc); ps['study']=ps.credible_set_uid.map(study)
def csw(d):
    return d.groupby('trait_scope').apply(lambda x: x.unresolved_context_mass.sum()/x.noncoding_pip_mass.sum()).round(5).to_dict()
print('credible-set-weighted unresolved, all:',csw(ps))
print('  EUR only:',csw(ps[ps.ancestry=='EUR']))
print('  no MVP:',csw(ps[~ps.study.str.startswith('MVP')]))
print('direct instances by study:',ps[ps.trait_scope=='tier1_direct_masld_pdff'].study.value_counts().to_dict())
