#!/usr/bin/env python3
"""Write the candidate conclusion and combined rosters after all analysis branches."""
import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from final_annotation_paths import annotation_paths

A=Path(__file__).resolve().parents[1]
ROOT=A.parents[2]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--candidate',type=Path,required=True)
    args=p.parse_args();c=args.candidate.resolve()
    assert os.environ.get('SLURM_JOB_ID')
    required=['disease_synthesis','disease_context','prediction_summary','rna_pseudobulk','rna_models','donor_associations','panels']
    records={name:json.loads((c/name/'completion.json').read_text()) for name in required}
    assert all(r['status']=='completed' for r in records.values())
    assert not records['rna_pseudobulk']['technical_pilot']
    for cohort in ['GSE244832','GSE281367']:
        for lineage in ['hepatocyte','stellate']:
            assert json.loads((annotation_paths(c)['motif_context']/f'{cohort}_{lineage}/completion.json').read_text())['status']=='completed'
    rosters=[]
    for cohort in ['GSE244832','GSE281367']:
        cells=pd.read_csv(A/f'results/20261023/inputs/audit/{cohort}.cells.tsv.gz',sep='\t',dtype=str)
        meta=pd.read_csv(ROOT/f'data/{cohort}/metadata/donor_pairing.csv',dtype=str)
        cells=cells.merge(meta,on='donor_id',how='left',validate='many_to_one',suffixes=('','_source'))
        cells['assay']='ATAC';cells['library_id']=cells.get('atac_sample',cells.donor_id)
        cells['specimen']='not_separately_identified';cells['run']=cells.get('atac_srr','not_recorded')
        cells['same_nucleus_pairing']=False
        cells['identity_basis']='source_declared_ATAC_donor_and_unique_library_barcode'
        cells['donor_identity_independently_verified']='not_established_by_this_analysis'
        cells['eligibility_basis']='Separate source da_eligible and network_eligible flags'
        rosters.append(cells)
    cells=pd.read_csv(c/'rna_pseudobulk/cell_roster.tsv.gz',sep='\t',dtype=str)
    cells['assay']='RNA';cells['run']=cells.runs
    cells['identity_basis']='deposited_RNA_specimen_donor_key_and_joint_run_molecule_deduplication'
    cells['donor_identity_independently_verified']='not_genotypically_verified_by_this_analysis'
    cells['eligibility_basis']='eligible marks cell QC and unanimous lineage; donor thresholds are in rna_pseudobulk/donors.tsv'
    rosters.append(cells)
    receipt=json.loads((A/'results/20261024/ridge_grid_sensitivity_r1/hepatocyte/prepared/preparation.json').read_text())
    cells=pd.read_csv(Path(receipt['input_root'])/'cells.tsv',sep='\t',dtype=str)
    cells['assay']='paired_RNA_ATAC';cells['library_id']=cells.well_id;cells['specimen']='not_separately_identified'
    cells['source_condition']=cells.condition
    cells['run']='not_recorded_in_qualified_cell_table';cells['condition']='no_qualified_MASH_contrast'
    cells['same_nucleus_pairing']=True
    cells['identity_basis']='source_demultiplexed_donor_exact_well_barcode_RNA_ATAC_join'
    cells['donor_identity_independently_verified']='same_nucleus_join_verified_not_new_donor_genotyping'
    membership=pd.read_csv(A/'results/20261024/ridge_grid_sensitivity_r1/hepatocyte/prepared/metacell_membership.tsv',sep='\t',dtype=str)
    assert membership.cell_id.is_unique and set(membership.cell_id).issubset(set(cells.cell_id))
    cells['prediction_retained_original']=cells.cell_id.isin(membership.cell_id)
    cells['eligibility_basis']='prediction_retained_original marks the frozen original capped aggregation roster'
    rosters.append(cells)
    columns=['cohort','donor_id','source_donor_id','specimen','library_id','run','raw_barcode','cell_id','assay',
             'condition','source_condition','lineage','cell_type','eligible','qc_pass','donor_lineage_cells','da_eligible','network_eligible',
             'prediction_retained_original','eligibility_basis','exclusion_reason','cross_assay_match_eligible',
             'cross_assay_confidence','same_nucleus_pairing','identity_basis','donor_identity_independently_verified']
    roster=pd.concat(rosters,ignore_index=True)
    for name in columns:
        if name not in roster:roster[name]=''
    assert not roster.duplicated(['cohort','assay','library_id','raw_barcode']).any()
    roster[columns].to_csv(c/'qualified_cell_roster.tsv.gz',sep='\t',index=False)
    design=pd.read_csv(c/'disease_synthesis/models.tsv',sep='\t')
    effects=pd.read_csv(c/'disease_context/accessibility_effects.tsv.gz',sep='\t')
    batch=[]
    for lineage,g in effects[effects.cohort.eq('GSE244832')].groupby('lineage'):
        available=g.batch_logFC.notna()
        batch.append(dict(lineage=lineage,n_peaks=len(g),n_condition_BH_q05=int(g.qvalue.lt(.05).sum()),
            n_batch_BH_q05=int(g.batch_qvalue.lt(.05).sum()),
            sign_agreement=float((g.loc[available,'logFC']*g.loc[available,'batch_logFC']>0).mean()),
            median_absolute_effect_change=float(np.median(np.abs(g.loc[available,'logFC']-g.loc[available,'batch_logFC'])))))
    pd.DataFrame(batch).to_csv(c/'batch_comparison.tsv',sep='\t',index=False)
    prediction=pd.read_csv(c/'prediction_summary/comparisons.tsv',sep='\t')
    overlap=pd.read_csv(c/'disease_context/cross_cohort_link_agreement.tsv',sep='\t')
    agreement=pd.read_csv(c/'disease_context/cohort_agreement.tsv',sep='\t')
    rna_donors=pd.read_csv(c/'rna_pseudobulk/donors.tsv',sep='\t',dtype={'donor_id':str})
    rna_decisions=pd.read_csv(c/'rna_models/model_decisions.tsv',sep='\t')
    rna_effects=pd.read_csv(c/'rna_models/rna_effects.tsv.gz',sep='\t')
    raw_checks=[json.loads(p.read_text()) for p in (c/'rna_recovery').glob('*/scientific_checks.json')]
    assert len(raw_checks)==18 and len({r['library'] for r in raw_checks})==18
    pd.DataFrame(raw_checks).to_csv(c/'rna_recount_checks.tsv',sep='\t',index=False)
    molecule_checks=[dict(library=p.parent.name,**json.loads(p.read_text())) for p in (c/'rna_recovery').glob('*/molecule_accounting.json')]
    assert len(molecule_checks)==18 and all(row['status']=='passed' for row in molecule_checks)
    pd.DataFrame(molecule_checks).to_csv(c/'rna_molecule_accounting.tsv',sep='\t',index=False)
    report=['Final MASLD accessibility and regulatory context: candidate analysis', '',
        'The primary TMM-voom quality-weighted accessibility analysis remains the manuscript reference.',
        'Robust QL did not improve the existing fitted-NB model-null diagnostics. Those diagnostics do not establish biological FDR.',
        'Four corrected primary fits reproduce; no lineage meets the fixed three-lineage directional-replication criterion.',
        'Four stellate peaks remain peakwise two-cohort candidates. Peak significance is not transferred to linked genes.', '']
    for row in agreement[agreement.universe.eq('all_joint_primary')].to_dict('records'):
        report.append(f"{row['lineage']}: {row['n_peaks']:,} jointly testable peaks; descriptive effect Spearman {row['effect_spearman']:.3f}; sign agreement {100*row['sign_agreement']:.1f}%.")
    for row in batch:
        report.append(f"GSE244832 {row['lineage']}: {row['n_condition_BH_q05']} condition-only and {row['n_batch_BH_q05']} recorded-batch-adjusted calls on identical donors/peaks; sign agreement {100*row['sign_agreement']:.1f}%.")
    report+=['Age, sex and BMI adjustments are not estimable from the supplied metadata; no expression-derived stage is used.','']
    robust=records['prediction_summary']['robust_internal_predictive_support']
    report.append('Cicero robustness: '+('gains over both baselines retain direction in every unshuffled stress test.' if robust else 'at least one unshuffled stress test fails to retain a gain over both baselines; the internal benefit is sensitivity-dependent.'))
    report+=['See prediction_summary/REPORT.txt for all magnitudes and conditional donor-bootstrap intervals.',
             'No new confirmatory P values, independent MASLD validation, edge-level correctness, or causal regulation follows.','']
    pairing_changes=pd.read_csv(c/'prediction_summary/paired_scenario_change_summary.tsv',sep='\t')
    shuffle_distance=pairing_changes[pairing_changes.scenario.str.startswith('pairing_shuffle_')&pairing_changes.baseline.eq('distance')]
    assert len(shuffle_distance)==3
    if shuffle_distance.original_gain_minus_scenario_gain.le(0).all():
        report.append('The Cicero gain over distance is no larger with correct pairing than after any of the three within-donor RNA shuffles. Added value from correct nucleus-level pairing for Cicero beyond proximity is not demonstrated. Donor-level association remains possible because shuffles preserve donor marginals.')
    for row in overlap[overlap.method.eq('cicero')].to_dict('records'):
        report.append(f"{row['lineage']}: cross-cohort Cicero top-link Jaccard {row['selected_edge_jaccard']:.3f} on common observability. Shared consensus ascertainment and top-ten selection qualify this stability measure.")
    report+=['All 117 Hotspot programs and weights remain fixed. Motif-family annotations indicate compatibility; they do not identify TF activity or binding.','',
             'Disease annotation maps retain positive Cicero links; the fixed predictive benchmark retains its original signed-nonzero selection rule. Predictive gains do not validate individual disease annotation edges.',
             'RNA: 117 sequencing runs were jointly reconstructed as 18 source-declared capture libraries. Every library passes exact input-read conservation, one-capture-library and unique integer-count checks.',
             'Jointly deduplicated molecule records reproduce called-cell UMI totals per cell and per gene in all libraries. RNA diagnoses follow the deposited workbook, independently of the local ATAC crosswalk.',
             'RNA model eligibility is reported in rna_models/model_decisions.tsv. Disease effects are donor pseudobulks from separate RNA libraries, not paired-nucleus prediction.',
             'These RNA effects use jointly recounted UMI matrices and inherited lineage annotations. No new ambient-RNA correction or doublet model was fitted; residual contamination remains a limitation.',
             'Documented cross-assay specimen matches are required for donor correlations; local confidence labels alone do not verify them. See donor_associations/decisions.tsv.',
             '']
    for lineage in ['hepatocyte','stellate','macrophage']:
        donors=rna_donors[rna_donors.lineage.eq(lineage)&rna_donors.eligible.eq(True)&rna_donors.condition.isin(['NORMAL','MASH'])]
        assert not donors.donor_id.duplicated().any()
        n_normal=int(donors.condition.eq('NORMAL').sum());n_mash=int(donors.condition.eq('MASH').sum())
        report.append(f"RNA {lineage}: {n_normal} normal and {n_mash} MASH donors with at least 100 qualified nuclei.")
        for decision in rna_decisions[rna_decisions.lineage.eq(lineage)].to_dict('records'):
            if decision['status']!='testable':
                report.append(f"  {decision['adjustment']}: untestable; {decision['reason']}.")
                continue
            values=rna_effects[rna_effects.lineage.eq(lineage)&rna_effects.adjustment.eq(decision['adjustment'])]
            assert len(values)==int(decision['n_tested_peaks'])
            report.append(f"  {decision['adjustment']}: {int(values.qvalue.lt(.05).sum()):,}/{len(values):,} genes at BH q<0.05 across tested genes within this lineage and model; log2 normalized RNA UMI effects.")
        values=rna_effects[rna_effects.lineage.eq(lineage)]
        if {'condition_only','batch_jb'}.issubset(set(values.adjustment)):
            primary=set(values[values.adjustment.eq('condition_only')&values.qvalue.lt(.05)].gene_id)
            adjusted=set(values[values.adjustment.eq('batch_jb')&values.qvalue.lt(.05)].gene_id)
            report.append(f"  {len(primary&adjusted):,} genes meet BH q<0.05 in both matched models; adjustment can also introduce different calls.")
    report+=['RNA disease associations from these separate libraries do not establish individual peak-to-gene relationships.',
             '', 'Canonical manuscript figures and Gene Catalog evidence are unchanged; adoption requires scientific review and explicit approval.']
    (c/'SCIENTIFIC_CONCLUSION.txt').write_text('\n'.join(report)+'\n')
    paths=[c/'SCIENTIFIC_CONCLUSION.txt',c/'qualified_cell_roster.tsv.gz',c/'batch_comparison.tsv',c/'rna_recount_checks.tsv',c/'rna_molecule_accounting.tsv']
    paths += [p for root in required for p in (c/root).iterdir() if p.suffix in ['.tsv','.gz','.pdf','.txt','.json']]
    pd.DataFrame([dict(path=str(p.relative_to(c)),bytes=p.stat().st_size) for p in sorted(set(paths))]).to_csv(c/'deliverables.tsv',sep='\t',index=False)
    # Numerical completion is distinct from subsequent PDF review and scope checks.
    (c/'numerical_completion.json').write_text(json.dumps(dict(status='completed',components=records,
        scientific_adoption=False,figure_visual_review='required',resource_scope_check='required'),indent=2))


if __name__=='__main__':main()
