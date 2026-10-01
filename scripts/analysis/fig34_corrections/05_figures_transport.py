"""Candidate panels, retrospective transport, and a comparison retaining native units."""
import sys
import numpy as np
import pandas as pd
from scipy import stats
from common import ROOT,HAC,ML,AXES,PROGRAM,FOCAL,SEED,read,write,fit,bh,plot_setup

def main(base):
    out=base/'figures_v4';out.mkdir(exist_ok=False);plt=plot_setup()
    md=read(base/'continuum/participant_histology.tsv')
    datasets=['GSE130970','GSE135251','GSE162694','GSE213621']
    sources=[('S4G',read(ML/'programs/hotspot/hotspot_participant_scores.tsv.gz').rename(columns={'program_uid':'feature_id'}),FOCAL),
             ('S4N',read(ML/'programs/nmf/nmf_participant_scores.tsv.gz').rename(columns={'outcome_id':'feature_id'}),None)]
    h=read(ML/'pathway_tf/pathways/hallmark/donor_scores.tsv.gz').rename(columns={'set_id':'feature_id','pathway_score':'outcome_z'})
    fixed=read(ML/'figures/source_tables/s3_six_hallmark_windows.tsv').set_id.unique()
    sources.append(('S3_hallmark_stage',h,{k:k.removeprefix('HALLMARK_').replace('_',' ') for k in fixed}))
    for prefix,ss,names in sources:
        if names is None:names={k:k for k in ss.feature_id.unique()}
        ss=ss[ss.feature_id.isin(names)].merge(md[['sample_id','fibrosis_group_reported','stage_resolution']],on='sample_id',validate='many_to_one')
        ss=ss[ss.dataset.isin(datasets)&ss.fibrosis_group_reported.ne('Not recorded')]
        tab=ss.groupby(['feature_id','dataset','fibrosis_group_reported']).outcome_z.agg(['mean','std','count']).reset_index()
        tab['biological_unit']='participant';tab['evidence_state']='descriptive';tab['interval']='pointwise_normal_approximation_mean'
        write(tab,out/f'{prefix}_source.tsv')
        for uid,name in names.items():
            fig,axs=plt.subplots(1,4,figsize=(7.2,1.95),sharey=True)
            for ax,co in zip(axs,datasets):
                order=['Control','F0–F1','F2','F3–F4'] if co=='GSE213621' else ['F0','F1','F2','F3','F4']
                z=tab[tab.feature_id.eq(uid)&tab.dataset.eq(co)].set_index('fibrosis_group_reported').reindex(order)
                ax.errorbar(range(len(order)),z['mean'],yerr=1.96*z['std']/np.sqrt(z['count']),fmt='o',ms=3,color='#C9265E',lw=.6)
                if co=='GSE213621':ax.plot(0,z['mean'].iloc[0],'o',ms=3,color='#9E9E9E');ax.set_facecolor('#F3F3F3')
                ax.set_title(co+('\nGrouped histology' if co=='GSE213621' else '\nExact stage'));ax.set_xticks(range(len(order)),order,rotation=25)
                ax.axhline(0,lw=.4,color='#9E9E9E')
            axs[0].set_ylabel('Mean stored score');fig.suptitle(name+' | Historical roster; identity conditional',fontsize=6);fig.tight_layout()
            fig.savefig(out/f'{prefix}_{uid}.pdf');plt.close(fig)
    # Main 4E retains the corrected stage contrasts, changing labels only.
    p=read(ROOT/'figures/candidates/fig4e-corrected-controls-20260925T175215Z/inputs/fig4e_bulk_tissue_state_transport.tsv')
    p['program_name']=p.program_uid.map(FOCAL);write(p,out/'fig4e_source.tsv')
    fig,ax=plt.subplots(figsize=(3.0,2.0))
    for uid,g in p.groupby('program_uid'):
        ax.plot(g.stage,g.effect,'o-',ms=3,lw=.7,label=FOCAL[uid],color='#C9265E' if FOCAL[uid]=='Stromal ECM' else '#2166AC')
    ax.axhline(0,lw=.5,color='#9E9E9E');ax.set_ylabel('Weighted member-gene log2FC');ax.set_xlabel('Cross-sectional contrast versus F0');ax.legend(frameon=False);ax.set_title('Historical roster; identity conditional');fig.tight_layout();fig.savefig(out/'fig4e_bulk_tissue_state_transport.pdf');plt.close(fig)
    # Preserve the actual main-panel maxT conjunction. Competitive tests are a
    # different hypothesis and remain separate source-table columns.
    src=ROOT/'figures/main/fig4_singlecell_programs/source_tables/current_candidate'
    a=read(src/'molecular_systems_atlas_20260824/overlays/community_stage_continuum_effects.tsv')
    labels=read(src/'molecular_systems_atlas_20260824/overlays/rendered_label_manifest.tsv')
    labels=labels[(labels.figure_id=='continuum_effect_map')&(labels.selection_rule=='15_largest_systems_fixed_before_outcome_overlay')]
    t=read(src/'molecular_systems_no_delta_20260825/source_tables/system_inference_strip.tsv')
    t=t[t.inference_column.eq('maxT')][['community_id','supported','criterion','metric_name','displayed_metric']]
    a=a.merge(t,on='community_id',validate='one_to_one');assert len(a)==43 and a.supported.sum()==18
    a['x_unit']='outcome_score_SD_per_cohort_SD_of_recorded_stage_code';a['y_unit']='outcome_score_SD_per_continuum_SD_adjusted_for_recorded_stage_category_and_sex'
    write(a,out/'fig4h_source.tsv')
    fig,ax=plt.subplots(figsize=(3.4,2.8))
    for sign,label,col in [(1,'Both positive','#C9265E'),(-1,'Both negative','#2166AC'),(0,'Opposite signs','#9E9E9E')]:
        mask=((a.stage_beta>0)&(a.continuum_beta>0)) if sign==1 else (((a.stage_beta<0)&(a.continuum_beta<0)) if sign==-1 else (a.stage_beta*a.continuum_beta<0))
        g=a[mask];ax.scatter(g.stage_beta,g.continuum_beta,c=col,s=15,label=label,edgecolors=np.where(g.supported,'black','none'),linewidths=.65)
    offsets={'S02':(-29,-14),'S04':(-26,13),'S07':(12,0),'S10':(-3,-23),'S13':(30,12),'S05':(0,13)}
    for r in a[a.community_id.isin(labels.community_id)].itertuples():
        offset=offsets.get(r.figure_label,(3,3))
        ax.annotate(r.figure_label,(r.stage_beta,r.continuum_beta),xytext=offset,textcoords='offset points',fontsize=5,arrowprops=dict(arrowstyle='-',lw=.35,color='#777777') if r.figure_label in offsets else None)
    ax.axhline(0,lw=.4,color='#9E9E9E');ax.axvline(0,lw=.4,color='#9E9E9E');ax.set_xlabel('Stage coefficient (score SD / cohort stage-code SD)');ax.set_ylabel('Stage-adjusted continuum coefficient\n(score SD / continuum SD)');ax.legend(frameon=False,loc='upper left',fontsize=5)
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([],[],marker='o',ls='',ms=3,color=col,label=label) for label,col in [('Both positive','#C9265E'),('Both negative','#2166AC'),('Opposite signs','#9E9E9E')]],frameon=False,loc='upper left',fontsize=5)
    ax.text(.98,.02,'Existing strict maxT outlines (18/43)',transform=ax.transAxes,ha='right',fontsize=5);ax.set_title('Historical roster; identity conditional');fig.tight_layout();fig.savefig(out/'fig4h_stage_continuum.pdf');plt.close(fig)
    # Cohort and omission results remain distinct from content-lineage follow-up.
    native=read(PROGRAM/'cohort_and_lodo_effects.tsv');native=native[native.program_uid.isin(FOCAL)].copy()
    native['unit']='native_Hotspot_score_per_ordinal_disease_grade';native['lineage_tested']='Hepatocytes'
    native['low']=native.beta-stats.t.ppf(.975,native.residual_df)*native.hc3_se;native['high']=native.beta+stats.t.ppf(.975,native.residual_df)*native.hc3_se
    write(native,out/'focal_native_cohort_and_omission.tsv')
    fig,axs=plt.subplots(1,2,figsize=(6,3.6))
    for ax,(uid,name) in zip(axs,FOCAL.items()):
        g=native[native.program_uid.eq(uid)&native.estimable].copy();g['label']=np.where(g.analysis_type.eq('cohort'),g.dataset_context,'Omit '+g.held_out_dataset.fillna(''))
        g=g[np.isfinite(g.low)&np.isfinite(g.high)]
        ax.errorbar(g.beta,np.arange(len(g)),xerr=[g.beta-g.low,g.high-g.beta],fmt='o',ms=2.5,lw=.6,color='#2166AC')
        ax.set_yticks(range(len(g)),[f'{r.label} (n={r.n_donors})' for r in g.itertuples()]);ax.axvline(0,color='#9E9E9E',lw=.5);ax.set_title(name);ax.set_xlabel('Native score / disease grade (HC3 95% CI)')
    fig.tight_layout();fig.savefig(out/'S4_native_replication_and_omission.pdf');plt.close(fig)
    # Previously evaluated cohorts: a declared retrospective exact-stage check.
    tr=read(HAC/'control_free_transfer/cft-v1-20260828T141625Z/per_sample_scores.tsv')
    rows=[]
    for co,g in tr.groupby('dataset'):
        if co not in ['GSE174478','GSE240729']:continue
        g=g[g.fibrosis_stage.notna()&g.fixed_projection.notna()].copy();assert not g.sample_id.duplicated().any()
        y=(g.fixed_projection-g.fixed_projection.mean())/g.fixed_projection.std(ddof=1)
        for adjustment in ['stage_only','stage_sex','stage_sex_age_complete_case']:
            keep=np.ones(len(g),bool);cols=[np.ones((len(g),1)),g[['fibrosis_stage']].to_numpy(float)]
            if adjustment!='stage_only':
                keep &= g.sex.notna();cols.append(pd.get_dummies(g.sex.fillna('missing'),drop_first=True).to_numpy(float))
            if adjustment=='stage_sex_age_complete_case':keep &=g.age.notna();cols.append(g[['age']].to_numpy(float))
            x=np.column_stack(cols);r=fit(y[keep],x[keep]);r.update(dataset=co,model=adjustment,unit='frozen_projection_cohort_SD_per_exact_stage',retrospective=True,assay='FFPE_bulk_RNA' if co=='GSE240729' else 'bulk_RNA');rows.append(r)
    trf=pd.DataFrame(rows)
    for _,idx in trf.groupby('model').groups.items():trf.loc[idx,'hc3_q_two_cohorts']=bh(trf.loc[idx,'hc3_p'],2)
    write(trf,out/'retrospective_exact_stage_transport.tsv')
    # Paired odds ratio uncertainty for the frozen expression-matched TREAT sensitivity.
    gs=ROOT/'figures/candidates/fig3g-coloc-r2-20260924T061535Z/source_tables'
    pair=read(gs/'fig3g_expression_matched_pairs.tsv');ga=pair.genetic_treat_supported;ma=pair.matched_control_treat_supported
    onlyg=int((ga&~ma).sum());onlym=int((~ga&ma).sum());both=int((ga&ma).sum());neither=int((~ga&~ma).sum())
    assert (len(pair),onlyg,onlym,both,neither)==(406,33,25,1,347)
    ci=stats.binomtest(onlyg,onlyg+onlym,.5).proportion_ci(method='exact')
    write(pd.DataFrame([dict(n_pairs=len(pair),both=both,neither=neither,genetic_only_TREAT=onlyg,matched_only_TREAT=onlym,paired_OR=onlyg/onlym,OR_low=ci.low/(1-ci.low),OR_high=ci.high/(1-ci.high),nominal_exact_McNemar_p=stats.binomtest(onlyg,onlyg+onlym,.5).pvalue,definition='expression_matched_TREAT_sensitivity_not_canonical_overlap')]),out/'fig3g_matched_uncertainty.tsv')
    sm=read(gs/'fig3g_matched_treat_summary.tsv').iloc[0];N=int(sm.n_jointly_testable);K=int(sm.n_genetically_anchored_joint);D=int(sm.n_canonical_deg_joint);O=int(sm.n_genetically_anchored_canonical)
    tab=np.array([[O,K-O],[D-O,N-K-D+O]])
    from scipy.stats.contingency import odds_ratio
    odds=odds_ratio(tab);oc=odds.confidence_interval()
    write(pd.DataFrame([dict(background=N,genetic=K,canonical=D,observed_overlap=O,expected_overlap=K*D/N,conditional_OR=odds.statistic,OR_low=oc.low,OR_high=oc.high,nominal_fisher_p=stats.fisher_exact(tab).pvalue,definition='canonical_jointly_testable_gene_background',evidence_state='conditional_on_historical_bulk_fit')]),out/'fig3g_canonical_overlap_uncertainty.tsv')
    scatter=read(gs/'fig3g_raw_treat_interface.tsv');assert len(scatter)==3622
    fig,ax=plt.subplots(figsize=(3.4,2.8))
    for state,col in [('Neither','#BBBBBB'),('Transcriptomics sig.','#C9265E'),('Genetics sig.','#2166AC'),('Both','#6A3D9A')]:
        g=scatter[scatter.state.eq(state)];ax.scatter(g.logFC,g.susie_pp4,s=3,c=col,alpha=.7,label=f'{state} ({len(g)})',linewidths=0)
    ax.set_xlabel('Disease-state log2FC');ax.set_ylabel('SuSiE PP.H4');ax.legend(frameon=False,fontsize=6);fig.tight_layout();fig.savefig(out/'fig3g_genetics_state_interface.pdf');plt.close(fig)
    (out/'COMPLETE').write_text('Candidate panels only; same sequence; preserved maxT hypothesis.\n')

if __name__=='__main__':
    from pathlib import Path
    main(Path(sys.argv[1]))
