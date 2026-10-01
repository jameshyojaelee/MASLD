"""Final candidate source tables and scientific-invariant checks; compute node."""
import sys,subprocess,hashlib,json,platform
import numpy as np
import pandas as pd
from common import ROOT,FOCAL,PROGRAM,read,write,plot_setup

def main(base):
    out=base/'source_tables_v2';out.mkdir(exist_ok=False);qa=base/'verification_v2';qa.mkdir(exist_ok=False)
    checks=[]
    def check(name,ok,detail):
        checks.append(dict(check=name,passed=bool(ok),detail=str(detail)))
        if not ok:write(pd.DataFrame(checks),qa/'failed_checks.tsv');raise AssertionError(name)
    for folder in ['continuum','protein','lineage_matched','coordination_v2','figures_v4','sensitivities','identity_sensitivity']:
        check(folder+'_complete',(base/folder/'COMPLETE').exists(),'producer completed')
    for file in ['continuum/inputs.tsv','lineage_matched/frozen_inputs.tsv']:
        for r in read(base/file).itertuples():
            p=__import__('pathlib').Path(r.path);sha=hashlib.sha256(p.read_bytes()).hexdigest();check('unchanged_'+p.name,sha==r.sha256,sha)
    member=read(PROGRAM/'program_membership_v2.tsv');check('117_programs_unchanged',member.program_uid.nunique()==117,member.program_uid.nunique())
    h=read(base/'continuum/participant_histology.tsv');group=h[h.dataset.eq('GSE213621')]
    check('combined_histology_no_exact_stage',group[group.fibrosis_group_reported.isin(['F0–F1','F3–F4','Control'])].fibrosis_stage_exact.isna().all(),'only documented F2 is exact')
    for f,col,limit in [('continuum/categorical_relabel_invariance.tsv','absolute_difference',1e-10),('continuum/stored_effect_reproduction.tsv','absolute_beta_difference',1e-8),('lineage_matched/raw_score_reproduction.tsv','max_abs_difference',1e-8),('coordination_v2/original_statistic_reproduction.tsv','max_abs_z_difference',1e-8)]:
        value=read(base/f)[col].max();check(f,value<limit,value)
    c=read(base/'continuum/continuum_cohort_effects.tsv')
    for (co,fam,axis,arm,comp),g in c.groupby(['dataset','family','axis','arm','composition']):
        check('family_'+str((co,fam,axis,arm,comp)),g.feature_id.nunique()==(117 if fam=='programs' else 43),len(g))
    l=read(base/'lineage_matched/lineage_scores.tsv.gz')
    keys=['donor','lineage','dataset','program_uid','minimum_cells','annotation_filter']
    rr=l[l.arm.eq('raw')&l.common_corrected_cells].merge(l[l.arm.eq('ambient_corrected_original_exposure')&l.common_corrected_cells],on=keys,suffixes=('_raw','_corrected'),validate='one_to_one')
    check('matched_raw_corrected_cells',rr.n_corrected_cells_raw.eq(rr.n_corrected_cells_corrected).all(),len(rr))
    pn=read(base/'protein/protein_pair_covariance.tsv');check('protein_same_46_cases',pn.n.eq(46).all(),'same stored participant frame for every model')
    # All exported inferential tables carry their interpretation with them.
    specs=[
      ('continuum/continuum_cohort_effects.tsv','continuum_cohort_effects.tsv','stored_outcome_SD_per_axis_SD','historical_source_library_roster_identity_unresolved','categorical_source_stage;inferred_sex;composition_when_flagged','117_programs_or_43_systems_per_cohort_axis_arm_model','conditional_historical_roster'),
      ('continuum/continuum_meta_effects.tsv','continuum_meta_effects.tsv','stored_outcome_SD_per_axis_SD','historical_source_library_roster_identity_unresolved','categorical_source_stage;inferred_sex;composition_when_flagged','117_programs_or_43_systems_per_axis_arm_model','conditional_historical_roster'),
      ('lineage_matched/lineage_cohort_effects.tsv','content_lineage_cohort_effects.tsv','projected_weighted_expression_score_difference','biological_donor','categorical_diagnosis_within_cohort','two_programs_per_contrast_arm_threshold_annotation','targeted_followup'),
      ('lineage_matched/lineage_meta_effects.tsv','content_lineage_meta_effects.tsv','projected_weighted_expression_score_difference','biological_donor','within_cohort_categorical_diagnosis_then_fixed_effect_meta','two_programs_primary_and_two_secondary_separate','targeted_followup'),
      ('protein/protein_pair_covariance.tsv','protein_pair_covariance.tsv','rank_residual_partial_correlation','participant','model_column_specifies_adjustment','three_pairs_per_model_and_twelve_pair_by_adjustment_tests','selection_conditioned'),
      ('coordination_v2/participant_bootstrap_summary.tsv','coordination_covariance.tsv','tanh_of_mean_pair_Fisher_z','source_participant_roster_identity_unresolved_not_edges','model_column_specifies_adjustment','pointwise_intervals_fixed_C3_not_multiplicity_corrected_tests','conditional_retrospective_covariance'),
      ('figures_v4/focal_native_cohort_and_omission.tsv','native_cohort_and_omission.tsv','native_Hotspot_score_per_ordinal_disease_grade','biological_donor','dataset_when_pooled;ordinal_diagnosis','original117_family_selection;cohort_omission_nominal_sensitivities','original_native_lineage'),
      ('figures_v4/fig4e_source.tsv','bulk_stage_transport.tsv','fixed_weighted_member_gene_log2FC','historical_source_library_roster_identity_unresolved','cohort;inferred_sex','gene_level23702_per_contrast;no_program_level_test','descriptive_transport'),
      ('sensitivities/captured_abundance_effects.tsv','captured_abundance.tsv','captured_lineage_log_odds','biological_donor','categorical_diagnosis_within_cohort','two_programs_per_contrast_arm_threshold_annotation','captured_fraction_not_tissue_fraction'),
      ('identity_sensitivity/conditional_continuum_effects.tsv','identity_conditional_continuum.tsv','stored_outcome_SD_per_axis_SD','one_library_per_provisional_identity_component','source_stage;inferred_sex;composition_when_flagged','117_programs_or_43_systems_per_cohort_axis_arm_model','conditional_identity_sensitivity'),
      ('sensitivities/age_complete_case_effects.tsv','age_sensitivity.tsv','stored_outcome_SD_per_complete_case_axis_SD','historical_source_library_roster_identity_unresolved','source_stage;inferred_sex;age_and_composition_when_flagged','117_programs_or_43_systems_per_axis_model','conditional_complete_case_sensitivity')]
    catalog=[]
    for source,name,unit,biounit,cov,family,state in specs:
        d=read(base/source);d['effect_unit']=unit;d['biological_unit']=biounit;d['covariate_definition']=cov;d['testing_family']=family;d['evidence_state']=np.where(d.estimable.eq(False),'untestable',state) if 'estimable'in d else state;d['source_file']=source
        write(d,out/name);catalog.append(dict(table=name,rows=len(d),effect_unit=unit,biological_unit=biounit,covariates=cov,testing_family=family,evidence_state=state))
    write(pd.DataFrame(catalog),out/'table_dictionary.tsv')
    # Render protein sensitivity with explicit selection conditioning, using
    # existing estimates; do not rerun the protein models for typography.
    plt=plot_setup();fig,axs=plt.subplots(1,3,figsize=(6.8,2.5),sharex=True)
    for ax,(pair,g) in zip(axs,pn.groupby('pair',sort=False)):
        ax.errorbar(g.rho,np.arange(len(g)),xerr=[g.rho-g.low,g.high-g.rho],fmt='o',ms=3,lw=.6,color='#2166AC');ax.set_yticks(range(len(g)),g.model.str.replace('_',' '));ax.set_title(pair);ax.axvline(0,color='#9E9E9E',lw=.5);ax.set_xlabel('Partial rank correlation (95% CI)')
    fig.suptitle('Selection-conditioned example; same 46 participants; pointwise bootstrap intervals',fontsize=6);fig.tight_layout();fig.savefig(base/'figures_v4/S3_protein_covariance.pdf');plt.close(fig)
    lf=read(base/'lineage_matched/lineage_cohort_effects.tsv')
    lf=lf[lf.minimum_cells.eq(50)&lf.annotation_filter.eq('all')&lf.contrast.str.startswith('primary')&lf.arm.ne('raw_full')]
    fig,axs=plt.subplots(1,2,figsize=(6.0,2.9))
    for ax,(uid,name) in zip(axs,FOCAL.items()):
        g=lf[lf.program_uid.eq(uid)]
        for j,r in enumerate(g.itertuples()):ax.errorbar(r.beta,j,xerr=[[r.beta-r.hc3_low],[r.hc3_high-r.beta]],fmt='o',ms=3,lw=.6,color='#9E9E9E' if r.arm=='raw' else '#C9265E')
        ax.set_yticks(range(len(g)),[f'{r.dataset} / '+('paired raw' if r.arm=='raw' else 'corrected')+f' (n={r.n})' for r in g.itertuples()]);ax.set_xlabel('Projected-score MASH − steatosis (HC3 95% CI)');ax.set_title(name);ax.axvline(0,color='#9E9E9E',lw=.5)
    fig.suptitle('Content-lineage follow-up; same correction-available cells and donors',fontsize=6);fig.tight_layout();fig.savefig(base/'figures_v4/S4_content_lineage_disease_effects.pdf');plt.close(fig)
    # Explicit missing category rows prevent a blank Control position from
    # being mistaken for a zero NMF loading.
    for stem in ['S4G','S4N','S3_hallmark_stage']:
        d=read(base/f'figures_v4/{stem}_source.tsv');rows=[]
        for (uid,co),g in d.groupby(['feature_id','dataset']):
            order=['Control','F0–F1','F2','F3–F4'] if co=='GSE213621' else ['F0','F1','F2','F3','F4']
            z=g.set_index('fibrosis_group_reported').reindex(order);z['feature_id']=uid;z['dataset']=co;z['count']=z['count'].fillna(0).astype(int);z['evidence_state']=np.where(z['count']>0,'descriptive_historical_roster','no_score_available');z['biological_unit']='historical_source_record_identity_unresolved';rows.append(z.reset_index())
        write(pd.concat(rows),out/f'{stem}_stage_display.tsv')
    pdfs=list((base/'figures_v4').glob('*.pdf'))+[base/'sensitivities/S4_focal_continuum_cohort_effects.pdf']
    proofs=qa/'pdf_text';proofs.mkdir();images=qa/'proofs';images.mkdir()
    manifest=[]
    for i,p in enumerate(pdfs):
        text=subprocess.check_output(['pdftotext','-raw',str(p),'-'],text=True)
        (proofs/(p.stem+'.txt')).write_text(text)
        check('PDF_text_'+p.name,len(text.strip())>20,len(text))
        if p.name.startswith(('S4G','S4N','S3_hallmark')):check('source_categories_'+p.name,all(s in text for s in ['Control','F0','F1','F2','F3','F4','Grouped histology']),'source grouped labels present')
        subprocess.run(['pdftoppm','-scale-to','1100','-singlefile','-png',str(p),str(images/p.stem)],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
        manifest.append(dict(panel=p.stem,path=str(p.resolve()),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),candidate_only=True))
    write(pd.DataFrame(manifest),qa/'figure_manifest.tsv');write(pd.DataFrame(checks),qa/'checks.tsv')
    (qa/'environment.txt').write_text(subprocess.check_output([sys.executable,'-m','pip','freeze'],text=True))
    (qa/'IDENTITY_LIMIT.txt').write_text('PASS of numerical invariants does not verify biological independence. Bulk historical-roster inference and GSE193066 evaluation remain conditional; final adoption is withheld.\n')
    (qa/'AUTOMATED_CHECKS_COMPLETE').write_text('Manual PDF review and scope check recorded separately.\n')

if __name__=='__main__':
    from pathlib import Path
    main(Path(sys.argv[1]))
