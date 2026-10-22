#!/usr/bin/env python3
"""Independent terminal table checks for the fixed H17 source-stage refit."""
import csv,gzip,hashlib,json,math,os
from pathlib import Path
ROOT=Path(os.environ['MASLD_PROJECT_ROOT']);OUT=Path(os.environ['H17_CHECK_ROOT']);OUT.mkdir(parents=True,exist_ok=False)
NEW=ROOT/'Analysis/Hypotheses/H17-source-stage-corrected-20261014_22237956';OLD=ROOT/'Analysis/Hypotheses/H17-20261009T015231Z/results_20261009T020140Z_job22127037';SRC=ROOT/'scripts/analysis/hypotheses_2026q4/H17_lipid_response_core'
def read(name,base=NEW):
 p=base/name;op=gzip.open if p.suffix=='.gz' else open
 with op(p,'rt') as h:return list(csv.DictReader(h,delimiter='\t'))
def quantile(values,p):
 a=sorted(values);i=(len(a)-1)*p;k=math.floor(i);return a[k]+(i-k)*(a[min(k+1,len(a)-1)]-a[k])
for name in ['STEP1_DONE','STEP2_DONE','STEP3_DONE']:assert (NEW/name).is_file(),name
assert (NEW/'input_sha256_before.tsv').read_bytes()==(NEW/'input_sha256_after.tsv').read_bytes()
for name in ['01_de_and_concordance.R','02_resampling.R','03_sensitivity_and_call.R']:assert (SRC/name).read_bytes()==(NEW/'scripts'/name).read_bytes()
insert='''
  qualified <- fread(Sys.getenv("H17_QUALIFIED_ROSTER"))
  stopifnot(!anyDuplicated(qualified$donor), all(donor_meta$sample %in% qualified$donor))
  j <- match(donor_meta$sample, qualified$donor)
  donor_meta[, disease_stage_coarse := qualified$disease_stage_coarse[j]]
  donor_meta[, exclude_stage_analysis := qualified$exclude_stage_analysis[j]]'''
lib=(NEW/'scripts/h17_lib.R').read_text();assert lib.count(insert)==1 and lib.replace(insert,'')==(SRC/'h17_lib.R').read_text()
assert 'H17_SEED <- 20260923L' in lib
counts={}
for name,key,n,col in [('S1_signperm_null.tsv.gz','perm',10000,'S1_null'),('S1_matchedset_null.tsv.gz','set',10000,'agreement'),('negctrl_invivo_label_permutation.tsv','perm',1000,'S1_agreement'),('negctrl_ipsc_label_swap.tsv','swap',100,'n_core'),('bootstrap_S2.tsv','boot',1000,'S2_spearman')]:
 r=read(name);assert len(r)==n and sorted(int(x[key]) for x in r)==list(range(1,n+1)),name
 assert all(math.isfinite(float(x[col])) for x in r),name;counts[name]=n
old={x['sample']:x for x in read('invivo_roster.tsv',OLD)};new={x['sample']:x for x in read('invivo_roster.tsv')};assert len(old)==59 and len(new)==55
assert set(old)-set(new)=={'SRR14602055','SRR14602056','SRR14602057','SRR14602058'} and not set(new)-set(old)
for key,x in new.items():
 for field in ['dataset','group','n_hep','n_runs','sex']:assert x[field]==old[key][field],(key,field)
assert sum(x['group']=='control' for x in new.values())==9 and len({x['dataset'] for x in new.values()})==2
assert len({x['participant'] for x in read('ipsc_sample_sheet.tsv')})==193
assert (NEW/'positive_control_ipsc.tsv').read_bytes()==(OLD/'positive_control_ipsc.tsv').read_bytes()
pc=read('positive_control_ipsc.tsv')
# Re-derive each assay-native BH family over all tested genes after filtering.
families={};de={}
for name in ['ipsc_de.tsv','invivo_de.tsv']:
 r=read(name);assert len({x['gene'] for x in r})==len(r);p=[float(x['P.Value']) for x in r];idx=sorted(range(len(p)),key=p.__getitem__);q=[0.]*len(p);run=1.
 for k in range(len(p)-1,-1,-1):
  j=idx[k];run=min(run,p[j]*len(p)/(k+1));q[j]=run
 assert max(abs(qi-float(x['q'])) for qi,x in zip(q,r))<1e-10,name
 families[name]=len(r);de[name]={x['gene']:x for x in r}
U=set(de['ipsc_de.tsv'])&set(de['invivo_de.tsv']);A=de['ipsc_de.tsv'];B=de['invivo_de.tsv'];sigA={g for g in U if float(A[g]['q'])<.05};sigB={g for g in U if float(B[g]['q'])<.05};same=sum(((float(A[g]['logFC'])>0)-(float(A[g]['logFC'])<0))==((float(B[g]['logFC'])>0)-(float(B[g]['logFC'])<0)) for g in sigA)
conc=read('concordance_primary.tsv')[0];assert len(U)==int(conc['n_U']) and len(sigA)==int(conc['n_sig_ipsc_in_U']) and len(sigB)==int(conc['n_sig_invivo_in_U']);assert abs(same/len(sigA)-float(conc['S1_agreement']))<1e-12
shared=len(sigA&sigB);expected=len(sigA)*len(sigB)/len(U);assert shared==int(conc['n_sig_both']) and abs(expected-float(conc['S3_expected_overlap']))<1e-10
core=sum(((float(A[g]['logFC'])>0)-(float(A[g]['logFC'])<0))==((float(B[g]['logFC'])>0)-(float(B[g]['logFC'])<0)) for g in sigA&sigB);assert core==int(conc['n_core'])
def ranks(values):
 idx=sorted(range(len(values)),key=values.__getitem__);out=[0.]*len(values);start=0
 while start<len(idx):
  end=start+1
  while end<len(idx) and values[idx[end]]==values[idx[start]]:end+=1
  rank=(start+1+end)/2
  for i in idx[start:end]:out[i]=rank
  start=end
 return out
ordered=sorted(U);ra=ranks([float(A[g]['logFC']) for g in ordered]);rb=ranks([float(B[g]['logFC']) for g in ordered]);mean=(len(U)+1)/2
rho=sum((a-mean)*(b-mean) for a,b in zip(ra,rb))/math.sqrt(sum((a-mean)**2 for a in ra)*sum((b-mean)**2 for b in rb));assert abs(rho-float(conc['S2_spearman']))<1e-12
choose=lambda n,k:math.lgamma(n+1)-math.lgamma(k+1)-math.lgamma(n-k+1)
upper=sum(math.exp(choose(len(sigA),k)+choose(len(U)-len(sigA),len(sigB)-k)-choose(len(U),len(sigB))) for k in range(shared,len(sigB)+1));assert abs(upper-float(conc['S3_hypergeom_p']))<1e-8

bs=[float(x['S2_spearman']) for x in read('bootstrap_S2.tsv')];lp=[float(x['S1_agreement']) for x in read('negctrl_invivo_label_permutation.tsv')];decision={x['criterion']:x for x in read('decision.tsv')}
# Re-derive controls and every prespecified support criterion rather than assume a verdict.
pc_up={}
for row in pc:
 g=row['gene'];expressed=g in A
 up=expressed and float(A[g]['q'])<.05 and float(A[g]['logFC'])>0
 assert (row['expressed']=='TRUE')==expressed and (row['up_sig']=='TRUE')==up,(g,'positive control')
 pc_up[g]=up
pc_pass=pc_up['PLIN2'] and sum(pc_up.values())>=math.ceil(sum(g in A for g in pc_up)/2)
assert all((x['overall_pass']=='TRUE')==pc_pass for x in pc)
sn=[float(x['S1_null']) for x in read('S1_signperm_null.tsv.gz')]
mn=[float(x['agreement']) for x in read('S1_matchedset_null.tsv.gz')]
agreement=same/len(sigA)
pupper=(1+sum(x>=agreement for x in sn))/(1+len(sn))
plower=(1+sum(x<=agreement for x in sn))/(1+len(sn))
matchedq95=quantile(mn,.95);labelq95=quantile(lp,.95)
for field,value in [('S1_p_upper',pupper),('S1_p_lower',plower),('S1_matchedset_q95',matchedq95)]:
 assert abs(float(conc[field])-value)<1e-12,(field,conc[field],value)
ci=[quantile(bs,.025),quantile(bs,.975)]
criteria={
 'positive_control_pass':pc_pass,
 'S1_agreement_gt_0.5_and_signperm_p_lt_0.05':agreement>.5 and pupper<.05,
 'S1_gt_matchedset_q95':agreement>matchedq95,
 'S1_gt_invivo_labelperm_q95':agreement>labelq95,
 'S2_bootstrap_CI_lower_gt_0':ci[0]>0,
 'core_nonempty_and_S3_p_lt_0.05':core>0 and upper<.05,
}
for criterion,met in criteria.items():
 assert decision[criterion]['met'] in ['TRUE','FALSE'] and (decision[criterion]['met']=='TRUE')==met,(criterion,decision[criterion],met)
verdict='untestable' if not pc_pass else 'yes' if all(criteria.values()) else 'no'
assert decision['CALL']['value']==verdict,('CALL',decision['CALL'],verdict)
assert (decision['S1_opposite_direction_lower_tail_p_lt_0.05']['met']=='TRUE')==(plower<.05)
assert f"[{ci[0]:.4f}, {ci[1]:.4f}]" in decision['S2_bootstrap_CI_lower_gt_0']['value']
support_count=sum(v for k,v in criteria.items() if k!='positive_control_pass')
swaps=read('negctrl_ipsc_label_swap.tsv');donor_null=read('negctrl_invivo_label_permutation.tsv')
negative_controls=dict(ipsc_swap_n_sig_median=quantile([int(x['n_sig_ipsc']) for x in swaps],.5),ipsc_swap_n_sig_max=max(int(x['n_sig_ipsc']) for x in swaps),ipsc_swap_core_max=max(int(x['n_core']) for x in swaps),donor_label_fraction_S1_le_observed=sum(float(x['S1_agreement'])<=agreement for x in donor_null)/len(donor_null),donor_label_fraction_rho_le_observed=sum(float(x['S2_spearman'])<=rho for x in donor_null)/len(donor_null),donor_label_fraction_core_ge_observed=sum(int(x['n_core'])>=core for x in donor_null)/len(donor_null))
report=dict(status='PASS',negative_controls=negative_controls,seed=20260923,participants_ipsc=193,in_vivo_donors=55,control_donors=9,case_donors=46,in_vivo_datasets=2,source_only_correction=True,resampling_counts=counts,BH_tested_gene_families=families,n_shared_genes=len(U),n_ipsc_significant_shared=len(sigA),same_sign_count=same,S1_agreement=same/len(sigA),donor_label_null_q95=quantile(lp,.95),S2_spearman=rho,S2_bootstrap_CI=ci,shared_significant=shared,expected_shared_significant=expected,nominal_overlap_p=float(conc['S3_hypergeom_p']),same_sign_core=core,positive_control=dict(prespecified_ipsc_pass=pc_pass,gene_up_significant=pc_up),independently_derived_criteria=criteria,S1_signperm_p_upper=pupper,S1_signperm_p_lower=plower,matchedset_q95=matchedq95,support_criteria_passed=support_count,terminal_verdict=verdict,sessionInfo_present=(NEW/'sessionInfo.txt').is_file())
assert report['sessionInfo_present'];(OUT/'verification.json').write_text(json.dumps(report,indent=2));(OUT/'CHECKED.txt').write_text('PASS: source-only annotation change; native seeds/models/resampling counts; biological units; BH gene families; every control/support criterion and terminal verdict independently rederived and matched to native decision.tsv.\n');print(json.dumps(report,indent=2))
