#!/usr/bin/env python3
"""Individual PDF panels for the final candidate, without manuscript writes."""
import argparse
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

GRAY='#9E9E9E'; MAGENTA='#C43C7C'; BLUE='#3E75AD'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':6,'axes.titlesize':6,'axes.labelsize':6,
 'xtick.labelsize':6,'ytick.labelsize':6,'legend.fontsize':6,'font.weight':'normal','axes.titleweight':'normal',
 'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'ps.fonttype':42})


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--candidate',type=Path,required=True)
    p.add_argument('--context',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args(); assert os.environ.get('SLURM_JOB_ID')
    args.output.mkdir(exist_ok=False)
    index=[]
    def save(fig,name,message,source):
        fig.tight_layout();fig.savefig(args.output/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
        index.append(dict(panel=name,message=message,source=str(source),file=name+'.pdf'))
    roster=pd.read_csv(args.candidate/'disease_synthesis/roster.tsv',sep='\t')
    # KEY MESSAGE: donor representation and nucleus depth differ by cohort and lineage.
    fig,ax=plt.subplots(figsize=(4.6,2.5));positions=[];labels=[]
    rng=np.random.default_rng(42)
    for i,((cohort,lineage),g) in enumerate(roster.groupby(['cohort','lineage'],sort=True)):
        for condition,color,offset in [('NORMAL',GRAY,-.15),('MASH',MAGENTA,.15)]:
            y=g.loc[g.condition.eq(condition),'n_cells']
            ax.scatter(i+offset+rng.uniform(-.06,.06,len(y)),y,s=9,c=color,label=condition if i==0 else None)
        positions.append(i);labels.append(cohort.replace('GSE','')+'\n'+lineage)
    ax.set_xticks(positions,labels);ax.set_yscale('log');ax.set_ylabel('ATAC nuclei per donor');ax.legend(frameon=False)
    save(fig,'01_donor_depth','Donors and nucleus depth define the inference units.',args.candidate/'disease_synthesis/roster.tsv')
    peaks=pd.read_csv(args.context/'cross_cohort_peak_evidence.tsv.gz',sep='\t')
    agreement=pd.read_csv(args.context/'cohort_agreement.tsv',sep='\t')
    for num,lineage in [(2,'hepatocyte'),(3,'stellate')]:
        # KEY MESSAGE: cohort agreement is descriptive and incomplete.
        g=peaks[peaks.lineage.eq(lineage)]
        fig,ax=plt.subplots(figsize=(3.2,2.8));ax.scatter(g.logFC_G244,g.logFC_G281,s=1,c=GRAY,alpha=.25,rasterized=True)
        called=g[g.peak_evidence_state.eq('two_cohort_peakwise_candidate')]
        if len(called):ax.scatter(called.logFC_G244,called.logFC_G281,s=12,c=MAGENTA)
        ax.axhline(0,c='black',lw=.4);ax.axvline(0,c='black',lw=.4)
        row=agreement[(agreement.lineage==lineage)&(agreement.universe=='all_joint_primary')].iloc[0]
        ax.set_title(f"{lineage}: n={len(g):,}; Spearman={row.effect_spearman:.3f}")
        ax.set_xlabel('GSE244832 MASH vs normal log2 fragment effect');ax.set_ylabel('GSE281367 log2 fragment effect')
        save(fig,f'{num:02d}_{lineage}_cohort_agreement','Peakwise agreement does not replace donor-level replication.',args.context/'cohort_agreement.tsv')
    effects=pd.read_csv(args.context/'accessibility_effects.tsv.gz',sep='\t')
    # KEY MESSAGE: recorded batch changes disease-effect estimates on identical donors and peaks.
    fig,ax=plt.subplots(figsize=(3.4,2.8))
    for lineage,color in [('hepatocyte',BLUE),('stellate',MAGENTA)]:
        g=effects[effects.cohort.eq('GSE244832')&effects.lineage.eq(lineage)&effects.batch_logFC.notna()]
        ax.scatter(g.logFC,g.batch_logFC,s=1,alpha=.2,c=color,rasterized=True,label=lineage)
    limits=ax.get_xlim();ax.plot(limits,limits,c=GRAY,lw=.6);ax.axhline(0,c='black',lw=.4);ax.axvline(0,c='black',lw=.4)
    ax.set_xlabel('Condition-only log2 fragment effect');ax.set_ylabel('Recorded-batch-adjusted effect');ax.legend(frameon=False,markerscale=4)
    save(fig,'04_batch_sensitivity','Batch adjustment is a sensitivity on matched donors and peaks.',args.context/'accessibility_effects.tsv.gz')
    # KEY MESSAGE: donor omissions expose effect sensitivity without changing biological n.
    influence=pd.read_csv(args.candidate/'disease_synthesis/influence.tsv',sep='\t')
    fig,ax=plt.subplots(figsize=(3.5,2.4))
    for i,((cohort,lineage),g) in enumerate(influence.groupby(['cohort','lineage'],sort=True)):
        ax.scatter(g.median_absolute_effect_change,np.repeat(i,len(g)),s=10,c=MAGENTA if lineage=='stellate' else BLUE)
    keys=list(influence.groupby(['cohort','lineage'],sort=True).groups)
    ax.set_yticks(range(len(keys)),[a.replace('GSE','')+' '+b for a,b in keys]);ax.set_xlabel('Median absolute change in log2 fragment effect')
    save(fig,'05_donor_influence','Each point is an independently refitted donor omission.',args.candidate/'disease_synthesis/influence.tsv')
    # KEY MESSAGE: selected Cicero links have limited cross-cohort overlap.
    overlap=pd.read_csv(args.context/'cross_cohort_link_agreement.tsv',sep='\t')
    fig,ax=plt.subplots(figsize=(3.4,2.2))
    for method,color,offset in [('distance',GRAY,-.08),('cicero',MAGENTA,.08)]:
        g=overlap[overlap.method.eq(method)].set_index('lineage').loc[['hepatocyte','stellate']]
        ax.scatter(g.selected_edge_jaccard,np.arange(2)+offset,s=16,c=color,label=method)
        for i,v in enumerate(g.selected_edge_jaccard):ax.annotate(f'{v:.3f}',(v,i+offset),xytext=(3,0),textcoords='offset points',va='center')
    ax.set_yticks([0,1],['hepatocyte','stellate']);ax.set_xlim(0,1.15);ax.set_xlabel('Top-link Jaccard on common observability');ax.legend(frameon=False,loc='lower center')
    save(fig,'06_link_agreement','Distance stability is expected from fixed geometry; overlap is not accuracy.',args.context/'cross_cohort_link_agreement.tsv')
    # KEY MESSAGE: program annotations retain the complete fixed family and original weights.
    programs=pd.read_csv(args.context/'fixed_program_link_coverage.tsv',sep='\t')
    fig,ax=plt.subplots(figsize=(3.4,2.4))
    for i,((cohort,lineage),g) in enumerate(programs[programs.method.eq('cicero')].groupby(['cohort','lineage'],sort=True)):
        ax.scatter(g.original_l1_weight_with_links,np.repeat(i,len(g))+rng.uniform(-.14,.14,len(g)),s=3,c=MAGENTA,alpha=.5)
    keys=list(programs[programs.method.eq('cicero')].groupby(['cohort','lineage'],sort=True).groups)
    ax.set_yticks(range(len(keys)),[a.replace('GSE','')+' '+b for a,b in keys]);ax.set_xlabel('Original fixed-program L1 weight with a Cicero link')
    save(fig,'07_fixed_program_coverage','All 117 programs are retained; annotation is not program significance.',args.context/'fixed_program_link_coverage.tsv')
    # KEY MESSAGE: example selection includes uncertain and batch-sensitive contexts.
    examples=pd.read_csv(args.context/'deterministic_examples.tsv',sep='\t')
    fig,ax=plt.subplots(figsize=(4.6,2.6))
    y=np.arange(len(examples));ax.errorbar(examples.logFC,y,xerr=[examples.logFC-examples.ci_lower,examples.ci_upper-examples.logFC],fmt='o',ms=3,c=MAGENTA,lw=.7)
    ax.axvline(0,c=GRAY,lw=.6);ax.set_yticks(y,[f'{r.display_category}: {r.gene_symbol}\n{r.peak_coordinate}' for r in examples.itertuples()]);ax.set_xlabel('MASH vs normal log2 fragment effect (95% moderated-t CI)')
    save(fig,'08_example_uncertainty','Deterministic examples preserve peak uncertainty and candidate status.',args.context/'deterministic_examples.tsv')
    pred=args.candidate/'prediction_summary/comparisons.tsv'
    if pred.exists():
        # KEY MESSAGE: modest gains must retain direction under fixed-fold stress tests.
        table=pd.read_csv(pred,sep='\t')
        fig,ax=plt.subplots(figsize=(4.4,3.4));names=list(dict.fromkeys(table.scenario))
        for comparator,color,offset in [('promoter_only',GRAY,-.1),('distance',MAGENTA,.1)]:
            g=table[table.baseline.eq(comparator)].set_index('scenario').loc[names]
            ax.errorbar(g.mean_gain,np.arange(len(names))+offset,xerr=[g.mean_gain-g.ci_low,g.ci_high-g.mean_gain],fmt='o',ms=3,c=color,lw=.6,label=comparator)
        ax.set_yticks(range(len(names)),names);ax.axvline(0,c='black',lw=.5);ax.set_xlabel('Baseline minus Cicero standardized MSE (95% donor interval)');ax.legend(frameon=False)
        save(fig,'09_prediction_robustness','Intervals condition on fitted predictions; internal development cohort.',pred)
        # KEY MESSAGE: shuffles test pairing while conserving donor RNA counts.
        changes=pd.read_csv(args.candidate/'prediction_summary/paired_scenario_change_summary.tsv',sep='\t')
        g=changes[changes.scenario.str.startswith('pairing_shuffle')&changes.baseline.eq('distance')]
        fig,ax=plt.subplots(figsize=(3.7,2.3));ax.errorbar(g.original_gain_minus_scenario_gain,range(len(g)),
            xerr=[g.original_gain_minus_scenario_gain-g.ci_low,g.ci_high-g.original_gain_minus_scenario_gain],fmt='o',ms=3,c=MAGENTA,lw=.7)
        ax.set_yticks(range(len(g)),g.scenario);ax.axvline(0,c=GRAY,lw=.6);ax.set_xlabel('Paired gain minus shuffled gain vs distance')
        save(fig,'10_pairing_controls','RNA shuffles preserve donor marginals and do not remove donor-level association.',args.candidate/'prediction_summary/paired_scenario_change_summary.tsv')
    rna=args.candidate/'rna_models/rna_effects.tsv.gz'
    if rna.exists() and rna.stat().st_size:
        # KEY MESSAGE: RNA corroboration retains source-library and donor eligibility limits.
        table=pd.read_csv(rna,sep='\t');table=table[table.adjustment.eq('condition_only')]
        fig,ax=plt.subplots(figsize=(3.6,2.4))
        for i,(lineage,g) in enumerate(table.groupby('lineage',sort=True)):
            ax.scatter(g.logFC,np.repeat(i,len(g))+rng.uniform(-.15,.15,len(g)),s=1,alpha=.15,c=MAGENTA,rasterized=True)
        names=list(table.groupby('lineage',sort=True).groups);ax.set_yticks(range(len(names)),names);ax.axvline(0,c=GRAY,lw=.6);ax.set_xlabel('Separate-library RNA MASH vs normal log2 UMI effect')
        save(fig,'11_RNA_disease_effects','RNA effects are donor-level separate-library corroboration.',rna)
    # KEY MESSAGE: robust QL did not improve the qualified fitted-NB calibration diagnostic.
    calibration=args.candidate.parents[1]/'20261023/da_diagnostic_review_qualified/null_calibration.tsv'
    null=pd.read_csv(calibration,sep='\t')
    null=null[null.run_arm.eq('full')&null.representation.eq('fragment')&null['filter'].eq('primary')]
    contexts=[('GSE244832','hepatocyte'),('GSE281367','hepatocyte'),('GSE244832','stellate'),('GSE281367','stellate')]
    fig,ax=plt.subplots(figsize=(3.8,2.5))
    for method,color,offset in [('voom_qw',GRAY,-.1),('edger_ql',MAGENTA,.1)]:
        g=null[null.method.eq(method)].set_index(['cohort','lineage']).loc[contexts]
        ax.errorbar(g.probability_any_BH_rejection,np.arange(len(contexts))+offset,
            xerr=[g.probability_any_BH_rejection-g.probability_any_BH_mc_lower,g.probability_any_BH_mc_upper-g.probability_any_BH_rejection],
            fmt='o',ms=3,c=color,lw=.7,label=method)
    ax.axvline(.05,c='black',lw=.5,ls=':');ax.set_yticks(range(4),[a.replace('GSE','')+' '+b for a,b in contexts]);ax.legend(frameon=False)
    ax.set_xlabel('Pr(any BH rejection), fitted-NB null; 1,000 draws')
    save(fig,'12_model_null_calibration','Model-null calibration is not biological FDR; binary NB nulls excluded.',calibration)
    pd.DataFrame(index).to_csv(args.output/'panels.tsv',sep='\t',index=False)
    (args.output/'completion.json').write_text(json.dumps(dict(status='completed',n_panels=len(index),individual_PDFs=True,scientific_adoption=False),indent=2))


if __name__=='__main__':main()
