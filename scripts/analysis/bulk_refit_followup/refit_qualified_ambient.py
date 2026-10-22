"""Refit existing frozen ambient scores after source-diagnosis qualification."""
import hashlib
import os
from pathlib import Path
import sys
import numpy as np
import pandas as pd

base = Path(os.environ['MASLD_PROJECT_ROOT'])
out = Path(os.environ['DIAGNOSIS_REFIT_ROOT'])
sys.path.insert(0, str(base/'Analysis/SingleCell/scripts/ambient_program_recalibration'))
from common import stage_fit, bh_adjust, evidence_state
old = base/'Analysis/SingleCell/candidates/ambient-program-recalibration-all117-candidate-2026-08-15-v2/results'
scores = pd.read_csv(old/'donor_program_scores.tsv.gz',sep='\t')
meta = pd.read_csv(out/'qualified_donor_metadata.tsv',sep='\t')
meta['exclude'] = meta.exclude_stage_analysis.astype(str).str.lower().isin(['true','t','1'])
# Preserve the original three-stage estimand. Source cirrhosis ordinal3 remains
# in metadata, but is not an eligible disease-stage contrast observation.
meta['stage_ordinal'] = meta.disease_stage_coarse.map({'Healthy':0.,'Steatosis':1.,'Steatohepatitis':2.})
effects = pd.read_csv(old/'ambient_program_effects.tsv',sep='\t')
reg = pd.read_csv(out/'hotspot/program_registry_v2.tsv',sep='\t').set_index('program_uid')
for i,row in effects.iterrows():
    frame = scores[scores.program_uid.eq(row.program_uid)].merge(meta,on='donor',validate='one_to_one')
    if row.analysis_universe == 'complete_case_common_universe':
        frame = frame[~frame.dataset.isin(['GSE189600','Liver_Atlas'])]
    frame['delta'] = frame.corrected-frame.raw
    for prefix in ['raw','corrected','delta']:
        try:
            result = stage_fit(frame,prefix)
            for key,value in result.items(): effects.loc[i,prefix+'_'+key] = value
            effects.loc[i,prefix+'_failure_reason'] = ''
        except RuntimeError as err:
            for col in effects.columns:
                if col.startswith(prefix+'_') and col!=prefix+'_failure_reason': effects.loc[i,col] = np.nan
            effects.loc[i,prefix+'_failure_reason'] = str(err)
    effects.loc[i,'n_model_donors'] = effects.loc[i,'raw_n_donors']
    effects.loc[i,'primary_selected_frozen'] = reg.loc[row.program_uid,'primary_selected']
    effects.loc[i,'robust_display_frozen'] = reg.loc[row.program_uid,'robust_display']
    effects.loc[i,'source_stability_frozen'] = reg.loc[row.program_uid,'source_stability']
for universe,idx in effects.groupby('analysis_universe').groups.items():
    assert len(idx)==117
    for prefix in ['raw','corrected','delta']:
        effects.loc[idx,prefix+'_qvalue'] = bh_adjust(effects.loc[idx,prefix+'_pvalue'],117)
        effects.loc[idx,prefix+'_hc3_qvalue'] = bh_adjust(effects.loc[idx,prefix+'_hc3_pvalue'],117)
effects['evidence_state']='sensitivity_only'
idx=effects.analysis_universe.eq('complete_case_common_universe')
effects.loc[idx,'evidence_state']=effects.loc[idx].apply(evidence_state,axis=1)
effects.to_csv(out/'ambient_program_effects.tsv',sep='\t',index=False)
meta.to_csv(out/'current_source_diagnosis_donor_roster.tsv',sep='\t',index=False)
inputs=[old/'donor_program_scores.tsv.gz',Path(os.environ['QUALIFIED_DONOR_METADATA'])]
inputs += [base/'Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot'/p for p in ['program_registry_v2.tsv','program_membership_v2.tsv','donor_program_scores_primary.tsv','donor_program_scores_equal_run.tsv']]
(out/'inputs.sha256').write_text(''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+str(p)+'\n' for p in inputs))
assert len(effects)==234
assert effects.loc[idx,'raw_n_donors'].max()==60
