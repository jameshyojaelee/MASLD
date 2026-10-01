"""Age complete-case sensitivity and captured-abundance uncertainty."""
import sys
import numpy as np
import pandas as pd
from scipy.linalg import helmert
from common import ROOT,HAC,ML,AXES,FOCAL,read,write,fit,bh,design,plot_setup

def main(base):
    out=base/'sensitivities';out.mkdir(exist_ok=False)
    md=read(base/'continuum/participant_histology.tsv')
    source=pd.read_csv(ROOT/'RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE162694/metadata/SraRunTable.csv')
    age=source[['Run','age']].rename(columns={'Run':'sample_id'})
    assert not age.sample_id.duplicated().any()
    md=md[md.dataset.eq('GSE162694')&~md.source_control&md.fibrosis_stage_exact.notna()].merge(age,on='sample_id',validate='one_to_one')
    md['age']=pd.to_numeric(md.age,errors='coerce');md=md.dropna(subset=['age'])
    a=read(AXES);md=md.merge(a[['sample_id','fixed_projection_raw','signature_pc1_raw']],on='sample_id',validate='one_to_one')
    f=pd.read_csv(ROOT/'Analysis/Deconvolution/results/GSE162694/GSE162694_bayesprism_proportions.tsv',sep='\t',index_col=0)
    ilr=pd.DataFrame(np.log(np.maximum(f.to_numpy(),1e-6))@(-helmert(16).T),index=f.index,columns=[f'ilr_{i:02d}' for i in range(15)])
    md=md.merge(ilr,left_on='sample_id',right_index=True,validate='one_to_one')
    p=read(ML/'programs/hotspot/hotspot_participant_scores.tsv.gz').rename(columns={'program_uid':'feature_id'});p['family']='programs'
    s=read(HAC/'molecular_systems_inference/hac-system-inference-20260824T232135Z/tables/participant_system_scores.tsv.gz').query("aggregation_method=='collection_balanced'").rename(columns={'community_id':'feature_id','system_score_z':'outcome_z'});s['family']='systems'
    dd=pd.concat([p[['sample_id','feature_id','outcome_z','family']],s[['sample_id','feature_id','outcome_z','family']]]).merge(md,on='sample_id',validate='many_to_one')
    rows=[]
    for axis in ['fixed_projection','signature_pc1']:
        for (fam,uid),g in dd.groupby(['family','feature_id']):
            g=g[np.isfinite(g.outcome_z)].copy();g['axis_z']=(g[axis+'_raw']-g[axis+'_raw'].mean())/g[axis+'_raw'].std(ddof=1)
            for comp in [False,True]:
                for adj in [False,True]:
                    x=design(g,composition=comp)
                    if adj:x=np.column_stack([x,g.age.to_numpy(float)-g.age.mean()])
                    r=fit(g.outcome_z,x);r.update(dataset='GSE162694',family=fam,feature_id=uid,axis=axis,composition=comp,age_adjusted=adj,unit='stored_outcome_SD_per_complete_case_axis_SD',covariates='categorical_exact_stage;inferred_sex'+(';15_composition_ILR' if comp else '')+(';source_age' if adj else ''));rows.append(r)
    rr=pd.DataFrame(rows)
    for (fam,axis,comp,adj),idx in rr.groupby(['family','axis','composition','age_adjusted']).groups.items():rr.loc[idx,'hc3_q']=bh(rr.loc[idx,'hc3_p'],117 if fam=='programs' else 43)
    write(rr,out/'age_complete_case_effects.tsv');write(md[['sample_id','age','fibrosis_stage_exact','inferred_sex']],out/'age_complete_case_participants.tsv')
    # Captured fractions use all six included atlas lineages, on exactly the
    # donors in each score contrast. This is not a tissue cell proportion.
    z=read(base/'lineage_matched/lineage_scores.tsv.gz');targets={list(FOCAL)[0]:'Fibroblasts',list(FOCAL)[1]:'Cholangiocytes'}
    rows=[]
    for (uid,co,minimum,ann,arm),g in z.groupby(['program_uid','dataset','minimum_cells','annotation_filter','arm']):
        g=g[g.lineage.eq(targets[uid]) & (g.common_corrected_cells | (arm=='raw_full'))]
        for label,ref in [('primary_steatohepatitis_vs_steatosis','Steatosis'),('secondary_steatohepatitis_vs_healthy','Healthy')]:
            d=g[g.disease.isin([ref,'Steatohepatitis'])];counts=d.disease.value_counts();x=np.column_stack([np.ones(len(d)),d.disease.eq('Steatohepatitis').to_numpy(float)])
            y=np.log((d.n_cells+.5)/(d.n_cells/d.captured_lineage_fraction-d.n_cells+.5))
            r=fit(y,x) if counts.get(ref,0)>=3 and counts.get('Steatohepatitis',0)>=3 else dict(n=len(d),estimable=False,beta=np.nan,hc3_p=np.nan,failure_reason='fewer_than_three_donors_per_group')
            r.update(program_uid=uid,dataset=co,minimum_cells=minimum,annotation_filter=ann,arm=arm,contrast=label,unit='captured_lineage_log_odds',n_reference=counts.get(ref,0),n_disease=counts.get('Steatohepatitis',0));rows.append(r)
    r=pd.DataFrame(rows)
    for _,idx in r.groupby(['dataset','minimum_cells','annotation_filter','arm','contrast']).groups.items():r.loc[idx,'hc3_q_two_programs']=bh(r.loc[idx,'hc3_p'],2)
    write(r,out/'captured_abundance_effects.tsv')
    coverage=z[z.arm.eq('raw_full')][['donor','lineage','dataset','disease','n_cells','n_corrected_cells','minimum_cells','annotation_filter','correction_cell_coverage','correction_complete']].drop_duplicates()
    coverage['correction_state']=np.where(coverage.n_corrected_cells==0,'no_corrected_cells',np.where(coverage.correction_complete,'complete','partial_cells_missing'))
    write(coverage,out/'correction_coverage.tsv')
    # Distinct assays and estimands remain in separate panels and tables.
    c=read(base/'continuum/focal_effects.tsv');c=c[c.axis.eq('fixed_projection')&c.arm.eq('disease_only')]
    plt=plot_setup();fig,axs=plt.subplots(1,2,figsize=(5.4,2.8))
    for ax,(uid,name) in zip(axs,FOCAL.items()):
        g=c[c.feature_id.eq(uid)]
        ax.errorbar(g.beta,np.arange(len(g)),xerr=[g.beta-g.hc3_low,g.hc3_high-g.beta],fmt='o',ms=3,lw=.6,color='#2166AC')
        ax.set_yticks(range(len(g)),[f'{r.dataset}, '+('composition adjusted' if r.composition else 'stage/sex')+f' (n={r.n})' for r in g.itertuples()]);ax.set_xlabel('Program SD / continuum SD (HC3 95% CI)');ax.set_title(name);ax.axvline(0,lw=.4,color='#9E9E9E')
    fig.tight_layout();fig.savefig(out/'S4_focal_continuum_cohort_effects.pdf');plt.close(fig)
    (out/'analysis_definition.txt').write_text('Age sensitivity: source-documented GSE162694 age, disease-only complete cases, matched baseline, frozen scores; full 117/43 BH families per axis and model. BMI not present in linked source table; GSE213621 source age/BMI unavailable, so no imputation or proxy adjustment. Captured abundance denominators are six captured lineages, not tissue fractions.\n')
    (out/'COMPLETE').write_text('Complete\n')

if __name__=='__main__':
    from pathlib import Path
    main(Path(sys.argv[1]))
