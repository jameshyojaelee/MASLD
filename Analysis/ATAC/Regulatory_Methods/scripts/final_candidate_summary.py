#!/usr/bin/env python3
"""Combine qualified disease estimates, all-donor links and fixed programs."""
import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from final_annotation_paths import annotation_paths

A=Path(__file__).resolve().parents[1]
ROOT=A.parents[2]
PROGRAM=ROOT/'Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot'
V3=ROOT/'Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-statistics-correction-20261022'


def coordinates(peaks):
    return peaks.chrom.astype(str)+':'+peaks.start0.astype(str)+'-'+peaks.end.astype(str)


def observable_peaks(directory):
    peaks=pd.read_csv(directory/'cicero/peaks.tsv',sep='\t')
    donors=pd.read_csv(directory/'cicero/donors.tsv',sep='\t',dtype=str)
    x=mmread(directory/'cicero/counts.mtx').T.tocsr()
    totals=np.asarray(x.sum(axis=1)).ravel()
    x=(sparse.diags(1e4/np.maximum(totals,1))@x).tocsr()
    positive=x.copy(); positive.data[:]=1
    logx=x.copy(); logx.data=np.log1p(logx.data)
    freq=donors.donor_id.value_counts()
    weights=np.array([1/freq[d] for d in donors.donor_id]); weights/=weights.sum()
    mean=np.asarray(weights@logx).ravel()
    square=logx.copy(); square.data**=2
    variance=np.asarray(weights@square).ravel()-mean**2
    keep=(np.asarray(weights@positive).ravel()>=.01)&(variance>1e-12)
    return set(coordinates(peaks.loc[keep]))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--candidate',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--allow-missing-rna',action='store_true')
    args=p.parse_args()
    assert os.environ.get('SLURM_JOB_ID')
    args.output.mkdir(exist_ok=False)
    disease=args.candidate/'disease_synthesis'
    effects=pd.read_csv(disease/'results.tsv.gz',sep='\t')
    stability=pd.read_csv(disease/'stability.tsv.gz',sep='\t')
    pk=['cohort','lineage','peak_coordinate']
    primary=effects[effects.adjustment.eq('condition_only')].copy()
    assert not primary.duplicated(pk).any()
    primary=primary.merge(stability,on=pk,how='left',validate='one_to_one')
    batch=effects[effects.adjustment.eq('batch_mm')][pk+['logFC','standard_error','pvalue','qvalue']].rename(columns={v:'batch_'+v for v in ['logFC','standard_error','pvalue','qvalue']})
    primary=primary.merge(batch,on=pk,how='left',validate='one_to_one')
    primary['batch_status']=np.where(primary.batch_logFC.notna(),'recorded_batch_adjustment_available','not_estimable_or_not_recorded')
    primary['batch_same_sign']=np.where(primary.batch_logFC.notna(),primary.logFC*primary.batch_logFC>0,np.nan)
    agreement=[]; peak_evidence=[]
    for lineage,rows in primary.groupby('lineage',sort=True):
        first=rows[rows.cohort.eq('GSE244832')].drop(columns='cohort')
        second=rows[rows.cohort.eq('GSE281367')].drop(columns='cohort')
        joined=first.merge(second,on=['lineage','peak_coordinate'],suffixes=('_G244','_G281'),validate='one_to_one')
        both=joined.qvalue_G244.lt(.05)&joined.qvalue_G281.lt(.05)
        same=joined.logFC_G244*joined.logFC_G281>0
        one=joined.qvalue_G244.lt(.05)^joined.qvalue_G281.lt(.05)
        joined['peak_evidence_state']=np.select([both&same,both&~same,one],['two_cohort_peakwise_candidate','discordant_peakwise','source_dependent'],default='indeterminate')
        joined['lineage_dynamic_main_criterion_met']=False
        for universe,frame in [('all_joint_primary',joined),('fixed_G281_discoveries',joined[joined.qvalue_G281.lt(.05)])]:
            agreement.append(dict(lineage=lineage,universe=universe,n_peaks=len(frame),
                effect_spearman=frame.logFC_G244.corr(frame.logFC_G281,method='spearman'),
                sign_agreement=float((frame.logFC_G244*frame.logFC_G281>0).mean()),
                interpretation='descriptive_peak_dependence_no_peak_independence_test'))
        peak_evidence.append(joined)
    pd.DataFrame(agreement).to_csv(args.output/'cohort_agreement.tsv',sep='\t',index=False)
    pe=pd.concat(peak_evidence,ignore_index=True)
    pe.to_csv(args.output/'cross_cohort_peak_evidence.tsv.gz',sep='\t',index=False)
    primary.to_csv(args.output/'accessibility_effects.tsv.gz',sep='\t',index=False)
    maps=[]; cover=[]; observability={}; shared_genes={}
    for cohort in ['GSE244832','GSE281367']:
        for lineage in ['hepatocyte','stellate']:
            directory=annotation_paths(args.candidate)['disease_links']/f'{cohort}_{lineage}'
            record=json.loads((directory/'completion.json').read_text()); assert record['status']=='completed'
            links=pd.read_csv(directory/'links.tsv.gz',sep='\t')
            cov=pd.read_csv(directory/'coverage.tsv',sep='\t')
            cov['cohort']=cohort; cov['lineage']=lineage; cover.append(cov)
            observability[(cohort,lineage)]=observable_peaks(directory)
            shared_genes[(cohort,lineage)]=set(cov.loc[cov.n_promoter_peaks.gt(0),'gene_id'])
            links=links.merge(primary,on=pk,how='left',validate='many_to_one',suffixes=('','_accessibility'))
            links['accessibility_testability']=np.where(links.logFC.notna(),'joint_primary_peak_tested','not_in_joint_primary_peak_universe')
            links['link_score_is_probability']=False
            maps.append(links)
    links=pd.concat(maps,ignore_index=True)
    links=links.merge(pe[['lineage','peak_coordinate','peak_evidence_state']],on=['lineage','peak_coordinate'],how='left',validate='many_to_one')
    links.peak_evidence_state=links.peak_evidence_state.fillna('not_jointly_tested')
    links['gene_regulatory_evidence']='accessibility_coupling_or_proximity_only'
    if (args.candidate/'rna_models/completion.json').exists():
        rna=pd.read_csv(args.candidate/'rna_models/rna_effects.tsv.gz',sep='\t')
        rna=rna[rna.adjustment.eq('condition_only')][['cohort','lineage','gene_id','logFC','standard_error','pvalue','qvalue']]
        rna=rna.rename(columns={v:'RNA_'+v for v in ['logFC','standard_error','pvalue','qvalue']})
        links=links.merge(rna,on=['cohort','lineage','gene_id'],how='left',validate='many_to_one')
        links['RNA_correspondence']=np.where(links.RNA_logFC.isna(),'untested_or_unavailable',
            np.where(links.logFC.isna(),'accessibility_untested',np.where(links.logFC*links.RNA_logFC>0,'same_direction','opposite_or_zero_direction')))
    else:
        assert args.allow_missing_rna,'RNA branch must be resolved before final summary'
        links['RNA_correspondence']='RNA_branch_unresolved_no_claim'
    links.to_csv(args.output/'peak_to_gene_context.tsv.gz',sep='\t',index=False)
    pd.concat(cover,ignore_index=True).to_csv(args.output/'gene_link_coverage.tsv.gz',sep='\t',index=False)
    overlap=[]
    for (lineage,method),group in links.groupby(['lineage','method']):
        a=group[group.cohort.eq('GSE244832')]; b=group[group.cohort.eq('GSE281367')]
        observable=observability[('GSE244832',lineage)]&observability[('GSE281367',lineage)]
        genes=shared_genes[('GSE244832',lineage)]&shared_genes[('GSE281367',lineage)]
        sets=[]
        for x in [a,b]:
            eligible=x[x.peak_coordinate.isin(observable)&x.gene_id.isin(genes)]
            sets.append(set(zip(eligible.gene_id,eligible.peak_coordinate)))
        common=sets[0]&sets[1]; union=sets[0]|sets[1]
        overlap.append(dict(lineage=lineage,method=method,n_common_observable_peaks=len(observable),
            n_G244_selected_observable=len(sets[0]),n_G281_selected_observable=len(sets[1]),n_shared_selected=len(common),
            selected_edge_jaccard=len(common)/len(union) if union else np.nan,
            interpretation='descriptive_top10_selection_agreement_on_common_observability_shared_consensus_ascertainment'))
    pd.DataFrame(overlap).to_csv(args.output/'cross_cohort_link_agreement.tsv',sep='\t',index=False)
    registry=pd.read_csv(PROGRAM/'program_registry_v2.tsv',sep='\t')
    membership=pd.read_csv(PROGRAM/'program_membership_v2.tsv',sep='\t')
    assert registry.program_uid.nunique()==117 and membership.program_uid.nunique()==117
    gene_axis=pd.read_csv(A/'results/20261023/preparation/genes.tsv',sep='\t')
    names=gene_axis.groupby('gene_name').gene_id.nunique()
    unique_symbols=set(names[names.eq(1)].index)
    membership['symbol_resolved_uniquely']=membership.mapped_symbol.isin(unique_symbols)
    membership.to_csv(args.output/'fixed_program_membership_annotation.tsv.gz',sep='\t',index=False)
    program_rows=[]
    for (cohort,lineage,method),group in links.groupby(['cohort','lineage','method']):
        linked=set(group.gene_symbol)&unique_symbols
        candidate=set(group.loc[group.peak_evidence_state.eq('two_cohort_peakwise_candidate'),'gene_symbol'])&unique_symbols
        for uid,mem in membership.groupby('program_uid',sort=True):
            supported=mem.mapped_symbol.isin(linked)&mem.symbol_resolved_uniquely
            eligible=mem.mapped_symbol.isin(candidate)&mem.symbol_resolved_uniquely
            program_rows.append(dict(program_uid=uid,cohort=cohort,lineage=lineage,method=method,
                n_fixed_genes=mem.mapped_symbol.nunique(),n_linked_genes=mem.loc[supported,'mapped_symbol'].nunique(),
                original_l1_weight_with_links=float(mem.loc[supported,'original_l1_weight'].sum()),
                n_genes_linked_to_peakwise_candidates=mem.loc[eligible,'mapped_symbol'].nunique(),
                no_program_redefinition=True,interpretation='fixed_program_annotation_not_program_significance'))
    pd.DataFrame(program_rows).to_csv(args.output/'fixed_program_link_coverage.tsv',sep='\t',index=False)
    # Display examples are selected from the complete table by a fixed ordering.
    candidates=links[links.logFC.notna()].copy()
    candidates['method_priority']=candidates.method.ne('cicero').astype(int)
    candidates['absolute_batch_change']=np.abs(candidates.batch_logFC-candidates.logFC)
    examples=[]
    for category,subset in [
        ('peakwise_candidate',candidates[candidates.peak_evidence_state.eq('two_cohort_peakwise_candidate')]),
        ('source_dependent',candidates[candidates.peak_evidence_state.eq('source_dependent')]),
        ('indeterminate',candidates[candidates.peak_evidence_state.eq('indeterminate')]),
        ('batch_sensitive',candidates[candidates.batch_logFC.notna()].sort_values('absolute_batch_change',ascending=False))]:
        if subset.empty: continue
        order=['method_priority','qvalue','lineage','peak_coordinate','gene_id','cohort']
        if category=='batch_sensitive':
            order=['absolute_batch_change']+order
            selected=subset.sort_values(order,ascending=[False]+[True]*(len(order)-1)).head(1).copy()
        else: selected=subset.sort_values(order).head(1).copy()
        selected['display_category']=category; examples.append(selected)
    pd.concat(examples,ignore_index=True).to_csv(args.output/'deterministic_examples.tsv',sep='\t',index=False)
    sources=[disease/'results.tsv.gz',disease/'checks.tsv',V3/'programs/program_atac_results.tsv',
             PROGRAM/'program_registry_v2.tsv',PROGRAM/'program_membership_v2.tsv']
    pd.DataFrame([dict(path=str(q),sha256=hashlib.sha256(q.read_bytes()).hexdigest()) for q in sources if q.exists()]).to_csv(args.output/'source_tables.tsv',sep='\t',index=False)
    (args.output/'completion.json').write_text(json.dumps(dict(status='completed',n_links=len(links),
        n_peakwise_candidates=int(pe.peak_evidence_state.eq('two_cohort_peakwise_candidate').sum()),
        n_programs=117,RNA_branch_resolved=(args.candidate/'rna_models/completion.json').exists(),
        scientific_adoption=False),indent=2))


if __name__=='__main__': main()
