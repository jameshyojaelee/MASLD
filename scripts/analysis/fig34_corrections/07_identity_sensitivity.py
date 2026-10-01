"""Conditional identity sensitivity; never adjudicates source-label conflicts.

Known same-person components with conflicting histology, diagnosis, cohort or
sex are excluded as a whole. A lexicographically first run is chosen within
concordant components, independently of target scores. Missing/possibly mixed
identity calls are excluded here, not silently declared independent people.
"""
import sys
import numpy as np
import pandas as pd
from scipy.linalg import helmert
from common import ROOT,HAC,ML,AXES,FOCAL,read,write,fit,bh,design,meta

def main(base):
    out=base/'identity_sensitivity';out.mkdir(exist_ok=False)
    ip=ROOT/'Analysis/MASLD_Model_Benchmark/executions/in-biopsy-allelic-a1-full-20260924T025858Z/identity'
    identity=read(ip/'sample_to_individual.tsv').rename(columns={'run':'sample_id'})
    pairs=read(ip/'kinship_related_pairs.tsv')
    m=read(base/'continuum/participant_histology.tsv').merge(identity,on='sample_id',how='left',validate='one_to_one')
    parent={s:s for s in m.sample_id}
    def root(s):
        parent.setdefault(s,s)
        if parent[s]!=s:parent[s]=root(parent[s])
        return parent[s]
    for a,b in pairs[pairs.relation.eq('same_individual')][['#IID1','IID2']].itertuples(index=False,name=None):
        aa,bb=root(a),root(b);parent[max(aa,bb)]=min(aa,bb)
    m['same_person_component']=m.sample_id.map(root)
    m['source_label_conflict']=False
    for _,g in m.groupby('same_person_component'):
        conflict=any(g[c].dropna().nunique()>1 for c in ['dataset','fibrosis_group_reported','group_binary','inferred_sex'])
        m.loc[g.index,'source_label_conflict']=conflict
    m['possibly_mixed']=m.expr_sex.eq('both')|m.other_allele_fraction_at_hom.gt(.02)
    m['reason']=np.select([m.source_label_conflict,m.possibly_mixed,m.individual_id.isna()|~m.passes_min_sites.fillna(False)],['conflicting_source_labels','possibly_mixed_identity_unresolved','insufficient_identity_evidence'],default='eligible')
    eligible=m[m.reason.eq('eligible')].sort_values('sample_id')
    chosen=set(eligible.drop_duplicates('same_person_component').sample_id)
    m.loc[m.reason.eq('eligible')&~m.sample_id.isin(chosen),'reason']='concordant_repeated_library'
    m['selected']=m.sample_id.isin(chosen);m['identity_state']='provisional_genotype_sensitivity_not_adopted'
    write(m,out/'roster_identity_disposition.tsv')
    census=m.groupby(['dataset','reason','source_control']).size().rename('n_libraries').reset_index();write(census,out/'identity_census.tsv')
    summary=m.groupby('dataset').agg(n_libraries=('sample_id','size'),same_person_components=('same_person_component','nunique'),selected_identity_qualified=('selected','sum'),conflicting_label_libraries=('source_label_conflict','sum'),possibly_mixed_libraries=('possibly_mixed','sum')).reset_index()
    summary['component_count_is_verified_donor_count']=False;write(summary,out/'roster_counts.tsv')
    scores=read(ML/'programs/hotspot/hotspot_participant_scores.tsv.gz').rename(columns={'program_uid':'feature_id'});scores['family']='programs'
    ss=read(HAC/'molecular_systems_inference/hac-system-inference-20260824T232135Z/tables/participant_system_scores.tsv.gz').query("aggregation_method=='collection_balanced'").rename(columns={'community_id':'feature_id','system_score_z':'outcome_z'});ss['family']='systems'
    d=pd.concat([scores[['sample_id','feature_id','outcome_z','family']],ss[['sample_id','feature_id','outcome_z','family']]]).merge(m[m.selected],on='sample_id',validate='many_to_one').merge(read(AXES)[['sample_id','fixed_projection_raw','signature_pc1_raw']],on='sample_id',validate='many_to_one')
    frac=[]
    for co in ['GSE162694','GSE213621']:
        f=pd.read_csv(ROOT/f'Analysis/Deconvolution/results/{co}/{co}_bayesprism_proportions.tsv',sep='\t',index_col=0)
        z=pd.DataFrame(np.log(np.maximum(f.to_numpy(),1e-6))@(-helmert(16).T),index=f.index,columns=[f'ilr_{i}' for i in range(15)]);z['sample_id']=f.index;frac.append(z)
    d=d.merge(pd.concat(frac),on='sample_id',how='left',validate='many_to_one');rows=[]
    for co in ['GSE162694','GSE213621']:
        c=d[d.dataset.eq(co)&d.fibrosis_stage.notna()]
        for arm in ['all_recorded_categories','disease_only']:
            z=c if arm=='all_recorded_categories' else c[~c.source_control]
            for axis in ['fixed_projection','signature_pc1']:
                for (fam,uid),g in z.groupby(['family','feature_id']):
                    g=g[np.isfinite(g.outcome_z)].copy();g['axis_z']=(g[axis+'_raw']-g[axis+'_raw'].mean())/g[axis+'_raw'].std(ddof=1)
                    for comp in [False,True]:
                        r=fit(g.outcome_z,design(g,composition=comp));r.update(dataset=co,family=fam,feature_id=uid,axis=axis,arm=arm,composition=comp,biological_unit='one_library_per_provisional_identity_qualified_component',evidence_state='conditional_sensitivity');rows.append(r)
    rr=pd.DataFrame(rows)
    for (co,fam,axis,arm,comp),idx in rr.groupby(['dataset','family','axis','arm','composition']).groups.items():rr.loc[idx,'hc3_q']=bh(rr.loc[idx,'hc3_p'],117 if fam=='programs' else 43)
    write(rr,out/'conditional_continuum_effects.tsv')
    # Mark identity limitations for the evaluation roster without applying
    # unreliable genotype calls to redefine the source participant crosswalk.
    ev=identity[identity.cohort.eq('GSE193066')].copy();ev['possibly_mixed']=ev.expr_sex.eq('both')|ev.other_allele_fraction_at_hom.gt(.02)
    write(ev,out/'GSE193066_identity_context.tsv')
    (out/'analysis_definition.txt').write_text('Conditional sensitivity only. No source-label adjudication or adopted cohort exclusion. Repeated libraries grouped by existing same-person pair graph; exclude whole components with discordant source histology/control/sex/cohort, missing identity calls, or possibly-mixed flags; choose lexicographic first eligible library per concordant component. Frozen scores are not rediscovered. GSE193066 identity uncertainty prevents interpreting its resampling as verified-donor replication.\n')
    (out/'COMPLETE').write_text('Identity limits documented; conditional sensitivity only.\n')

if __name__=='__main__':
    from pathlib import Path
    main(Path(sys.argv[1]))
