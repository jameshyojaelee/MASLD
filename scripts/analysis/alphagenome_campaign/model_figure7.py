#!/usr/bin/env python3
"""Individual candidate Figure 7 panels from completed development comparisons.

# KEY MESSAGE: These initial reporter and donor-difference comparisons do not
# establish an improvement over zero molecular differences on held examples.
"""
import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd

FONT=6
GRAY="#9E9E9E"
BLUE="#1565C0"
MAGENTA="#C9265E"


def setup():
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":FONT,"axes.titlesize":FONT,
        "axes.labelsize":FONT,"xtick.labelsize":FONT,"ytick.labelsize":FONT,"legend.fontsize":FONT,
        "axes.titleweight":"normal","font.weight":"normal","pdf.fonttype":42,"ps.fonttype":42,
        "axes.spines.top":False,"axes.spines.right":False,"axes.linewidth":0.5,
        "xtick.major.width":0.5,"ytick.major.width":0.5,"legend.frameon":False})


def common(ax,title,xlabel):
    ax.axvline(0,color=GRAY,lw=0.7,ls="--",zorder=0)
    ax.set_title(title,loc="left",pad=7)
    ax.set_xlabel(xlabel)
    ax.tick_params(length=2,pad=2)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=4))


def main(args):
    args.out.mkdir(parents=True,exist_ok=False)
    setup()
    values=[]
    paths=[args.lx2/"performance.tsv",args.lx2/"paired_intervals.tsv",args.lx2/"population_coverage.tsv",
           args.donor/"metrics.tsv",args.donor/"paired_uncertainty.tsv",args.donor/"receipt.json"]
    performance=pd.read_csv(paths[0],sep="\t")
    intervals=pd.read_csv(paths[1],sep="\t")
    population=pd.read_csv(paths[2],sep="\t")
    arms=[("direct_interaction","Direct interaction"),("separate_conditions","Separate condition heads"),
          ("shared_condition_residual","Shared + condition residual")]
    for cell,label in (("LX2","LX-2"),("HepG2","HepG2")):
        panel=f"Figure7_{cell}_reporter_interactions"
        fig,ax=plt.subplots(figsize=(3.9,2.1))
        base=float(performance.loc[(performance.cell_line == cell)&(performance.arm == "zero_interaction"),"RMSE_log2_activity_interaction"].iloc[0])**2
        for index,(arm,text) in enumerate(arms):
            row=performance.loc[(performance.cell_line == cell)&(performance.arm == arm)].iloc[0]
            interval=intervals.loc[(intervals.cell_line == cell)&(intervals.arm == arm)].iloc[0]
            point=base-float(row.RMSE_log2_activity_interaction)**2
            lo,hi=interval.MSE_improvement_vs_zero_low95,interval.MSE_improvement_vs_zero_high95
            ax.plot([lo,hi],[index,index],color=BLUE,lw=0.7)
            ax.plot(point,index,"o",color=BLUE,ms=3)
            values.append({"panel":panel,"row":text,"context":cell,"point":point,"low95":lo,"high95":hi,
                           "n_constructs":int(row.constructs),"n_blocks":int(row.locus_blocks),"experimental_replicates":4})
        ax.set_yticks(range(3),[x[1] for x in arms]);ax.invert_yaxis()
        ax.set_ylim(2.5,-0.5)
        common(ax,f"{label}: treatment interaction","MSE reduction vs zero interaction (log2 activity)²")
        pop=population.loc[population.cell_line == cell].iloc[0]
        note=f"{int(pop.matched_four_pair_constructs):,} of {int(pop.full_source_constructs):,} constructs; 239 locus blocks\n4 paired experimental replicates; development comparison"
        fig.text(0.03,0.04,note,fontsize=FONT)
        fig.subplots_adjust(left=0.38,right=0.98,top=0.82,bottom=0.35)
        fig.savefig(args.out/(panel+".pdf"));plt.close(fig)
    donor_metrics=pd.read_csv(paths[3],sep="\t")
    donor_ci=pd.read_csv(paths[4],sep="\t")
    donor_ci=donor_ci.loc[donor_ci.endpoint == "donor_difference_mse_log2cpm"].copy()
    tasks=[("held_donors_trained_regions","Held donors, trained regions",BLUE,"o",-0.12),
           ("held_donors_held_regions","Held donors and regions",MAGENTA,"s",0.12)]
    recipes=[("all_rna:shared_profile","Profile objective"),("all_rna:shared_training_mean_residual","Training-mean residual"),
             ("all_rna:shared_pairwise","Donor-pair objective"),("all_rna:additive","Additive sequence + RNA"),
             ("all_rna:rna_only_region_specific","Region-specific RNA")]
    fig,ax=plt.subplots(figsize=(4.4,2.7))
    panel="Figure7_donor_objectives"
    for evaluation,label,color,marker,offset in tasks:
        for index,(model,text) in enumerate(recipes):
            selected=donor_ci.loc[(donor_ci.evaluation == evaluation)&(donor_ci.model == model)]
            if selected.empty:
                continue
            row=selected.iloc[0]
            point,lo,hi=-row.mse_difference,-row.ci95_high,-row.ci95_low
            ax.plot([lo,hi],[index+offset,index+offset],color=color,lw=0.7)
            ax.plot(point,index+offset,marker,color=color,ms=3,label=label if index == 0 else None)
            values.append({"panel":panel,"row":text,"context":evaluation,"point":point,"low95":lo,"high95":hi,
                           "n_donors":99,"n_regions":128 if evaluation == tasks[0][0] else 64})
    ax.set_yticks(range(len(recipes)),[x[1] for x in recipes]);ax.invert_yaxis()
    common(ax,"Donor-difference prediction: bounded linear pilot","MSE reduction vs zero donor difference\n(log2[1 + H3K27ac CPM])²")
    fig.legend(*ax.get_legend_handles_labels(),loc="upper left",bbox_to_anchor=(0.35,0.98),borderaxespad=0,handlelength=1.2)
    fig.text(0.03,0.025,"99 participants; 128 regions (64 held regions)\nRegion-specific RNA has no matched held-region comparison.",fontsize=FONT)
    fig.subplots_adjust(left=0.35,right=0.98,top=0.75,bottom=0.32)
    fig.savefig(args.out/(panel+".pdf"));plt.close(fig)
    exclusions=[("all_rna","All eligible RNA"),("overlapping_genes_excluded","Overlapping genes excluded"),
        ("within_1mb_excluded","RNA within ±1 Mb excluded"),("chromosome_excluded","Target chromosome excluded"),
        ("random_overlap","Matched random: overlap"),("random_1mb","Matched random: ±1 Mb"),
        ("random_chromosome","Matched random: chromosome"),("technical_only","Technical inputs only"),
        ("context_shuffled","RNA context shuffled")]
    panel="Figure7_donor_context_exclusions"
    fig,ax=plt.subplots(figsize=(4.4,3.4))
    for evaluation,label,color,marker,offset in tasks:
        for index,(key,text) in enumerate(exclusions):
            selected=donor_ci.loc[(donor_ci.evaluation == evaluation)&(donor_ci.model == key+":shared_pairwise")]
            if len(selected) != 1:
                raise ValueError(f"Missing/duplicate donor comparison {evaluation}/{key}")
            row=selected.iloc[0];point,lo,hi=-row.mse_difference,-row.ci95_high,-row.ci95_low
            ax.plot([lo,hi],[index+offset,index+offset],color=color,lw=0.7)
            ax.plot(point,index+offset,marker,color=color,ms=3,label=label if index == 0 else None)
            values.append({"panel":panel,"row":text,"context":evaluation,"point":point,"low95":lo,"high95":hi,
                           "n_donors":99,"n_regions":128 if evaluation == tasks[0][0] else 64})
    ax.set_yticks(range(len(exclusions)),[x[1] for x in exclusions]);ax.invert_yaxis()
    common(ax,"Donor-pair objective: RNA exclusion sensitivities","MSE reduction vs zero donor difference\n(log2[1 + H3K27ac CPM])²")
    fig.legend(*ax.get_legend_handles_labels(),loc="upper left",bbox_to_anchor=(0.42,0.98),borderaxespad=0,handlelength=1.2)
    fig.text(0.03,0.02,"99 participants; simple linear pilot; descriptive 95% intervals\nOverlapping genes are a sensitivity, not validated target removal.",fontsize=FONT)
    fig.subplots_adjust(left=0.42,right=0.98,top=0.79,bottom=0.25)
    fig.savefig(args.out/(panel+".pdf"));plt.close(fig)
    pd.DataFrame(values).to_csv(args.out/"plotted_values.tsv",sep="\t",index=False)
    manifest={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    (args.out/"sources.json").write_text(json.dumps(manifest,indent=2)+"\n")
    reporter_populations="; ".join(
        f"{r.cell_line}: {int(r.matched_four_pair_constructs):,} of {int(r.full_source_constructs):,} source constructs"
        for r in population.itertuples())
    (args.out/"CAPTIONS.md").write_text("""# Candidate Figure 7 panels: initial development comparisons

These panels are candidates and are not adopted Resource figures. Positive
values mean lower mean squared error than the zero-difference comparator.
Intervals are descriptive 95% intervals with no confirmatory significance claim.

**Reporter interactions.** Frozen local AlphaGenome 2,048-bp variant-bin
representations with forward/reverse-complement averaging; REPORTER_POPULATIONS, 239 held locus
blocks, and four paired experimental replicates per condition. All feature
compression, scaling and penalty selection are training-only within nested
development folds. Intervals use 2,000 paired locus and jointly resampled
experimental-replicate draws (seed 20260915), retaining shared-control
covariance. Fitting/selection uncertainty is not included. LX-2 TGFβ and HepG2
palmitate/oleate interaction targets are cell-line reporter measurements; they
do not establish participant fibrosis or treatment response. None of the
three fitted arms has lower point-estimate error than zero interaction.

**Donor objectives and exclusions.** This is a 99-participant, 128-region
GSE267145 pilot with simple linear RNA and counted-interval dinucleotide/width
representations. The crossed evaluation holds 64 output regions. All output
errors use log2(1 + H3K27ac CPM), without evaluation-fitted scaling. Intervals
use 1,000 paired donor-by-tile bootstrap draws. Region-specific RNA is shown
only at trained output regions. The exclusion ladder tests overlapping genes,
nearby RNA and chromosome-local RNA, matched random exclusions, technical
inputs and shuffled context. Validated target links and batch metadata were
unavailable, so these are not validated target-removal or batch-stratified
controls. Counted-interval sequence coordinates carry a 1-bp uncertainty.
The pilot does not compare full nonlinear donor-context models or establish
trans regulation. Protected outcomes remain unopened.
""".replace("REPORTER_POPULATIONS",reporter_populations))
    print(f"Wrote four candidate PDFs and {len(values)} plotted-value rows to {args.out}")


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lx2",type=Path,required=True)
    parser.add_argument("--donor",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    main(parser.parse_args())
