#!/usr/bin/env python3
"""Candidate adaptation panel from independently reconstructed development errors.

# KEY MESSAGE: Six adapter gains survive paired uncertainty against their fixed
# controls in one development split, while native-model superiority is unresolved.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd


def label(arm):
    if arm.startswith("partial"):
        return "Partial q/v, " + ("3e-6" if arm.endswith("3e-06") else "1e-5")
    _, rank, blocks, pool = arm.split("_")
    return f"Rank {rank[1:]}, last {blocks[4:]}, {pool}"


def main(args):
    if not os.environ.get("SLURM_JOB_ID"):
        raise SystemExit("Generate candidate figures on a compute node")
    args.out.mkdir(parents=True, exist_ok=False)
    check = json.loads(args.check.read_text())
    assert check["status"] == "pass" and check["contrasts"] == 20
    source = pd.read_csv(args.source,sep="\t")
    chromosome = source.loc[source.uncertainty_unit.eq("whole_chromosomes")].sort_values("RMSE_improvement_beta",ascending=False)
    assert len(chromosome) == 10 and chromosome.units.eq(5).all()
    bins = source.loc[source.uncertainty_unit.eq("historical_1mb_bins")].set_index("arm").loc[chromosome.arm]
    np.testing.assert_allclose(chromosome.RMSE_improvement_beta,bins.RMSE_improvement_beta,rtol=0,atol=1e-14)
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":6,"axes.titlesize":6,
        "axes.labelsize":6,"xtick.labelsize":6,"ytick.labelsize":6,"legend.fontsize":6,
        "axes.titleweight":"normal","font.weight":"normal","pdf.fonttype":42,"ps.fonttype":42,
        "axes.spines.top":False,"axes.spines.right":False,"axes.linewidth":.5,
        "xtick.major.width":.5,"ytick.major.width":.5,"legend.frameon":False})
    fig, ax = plt.subplots(figsize=(4.9,3.65))
    ax.axvline(0,color="#9E9E9E",lw=.7,ls="--")
    values = []
    for index, row in enumerate(chromosome.itertuples(index=False)):
        previous = bins.loc[row.arm]
        for unit, record, offset, marker, color in (
            ("524 historical 1-Mb bins",previous,-.13,"o","#1565C0"),
            ("5 chromosomes",row,.13,"s","#C9265E")):
            low, high = record.RMSE_improvement_beta_low95,record.RMSE_improvement_beta_high95
            point = record.RMSE_improvement_beta
            ax.plot([low,high],[index+offset,index+offset],color=color,lw=.8)
            ax.plot(point,index+offset,marker,color=color,ms=3,label=unit if index == 0 else None)
            values.append(dict(arm=row.arm,label=label(row.arm),uncertainty_unit=unit,
                point=point,low95=low,high95=high,source_units="source_FastQTL_ALT_dosage_beta",
                exploratory_BH_q=record.exploratory_BH_q, rows=6834))
        ax.text(1.025,index,f"{row.RMSE_improvement_beta:.4f}",transform=ax.get_yaxis_transform(),va="center")
    ax.set_yticks(range(10),[label(arm) for arm in chromosome.arm])
    ax.set_ylim(9.6,-.6)
    ax.set_xlabel("RMSE reduction vs matching fixed control\n(source accessibility β units)")
    ax.set_title("Adaptation at fixed training exposure",loc="left",pad=7)
    ax.tick_params(length=2,pad=2)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.text(1.025,1.015,"Reduction",transform=ax.transAxes,va="bottom")
    fig.legend(*ax.get_legend_handles_labels(),loc="upper left",bbox_to_anchor=(.37,.995),
        handlelength=1.2,borderaxespad=0)
    fig.text(.025,.025,"6,834 variants; one development split and seed; 5,000 training steps\n"
        "Conditional 95% intervals; full native-predictor comparison pending.")
    fig.subplots_adjust(left=.37,right=.90,top=.80,bottom=.22)
    fig.savefig(args.out / "Figure7_adaptation_development.pdf")
    plt.close(fig)
    pd.DataFrame(values).to_csv(args.out/"plotted_values.tsv",sep="\t",index=False)
    (args.out/"sources.json").write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest()
        for p in (args.source,args.check)},indent=2)+"\n")
    (args.out/"CAPTIONS.md").write_text("""# Candidate Figure 7: adaptation development comparison

All eight adapter and two partial query/value update recipes are shown against
their pooling-matched fixed 2,048-bp shared-head controls. Positive values mean
lower RMSE in the source FastQTL accessibility slope per ALT dosage. Rows are
ordered by observed RMSE reduction; this does not select a final recipe.

The comparison uses the same 6,834 variants in development fold 1, 5,000 training
steps, microbatch 4, and seed 1103. Training uses folds 2–4; fold 0 is not used in
these adaptation evaluations. Intervals use 2,000 paired bootstrap draws over
524 historical 1-Mb bins (seed 20260915) or five chromosomes (seed 20260916).
The historical bins are not verified independent LD blocks, and sequence
windows overlap across bins. Five chromosome units provide a limited
sensitivity analysis. Variants and chromosomes are not biological participants.

This precision analysis was specified after aggregate pilot inspection.
Native-unit MSE contrasts form a fixed exploratory family of ten, adjusted
separately with BH for each resampling unit. Six of eight adapter MSE-gain
intervals exclude zero and have BH q < 0.05 under both units. No asterisks or
confirmatory claim are attached to this post-pilot comparison. Neither partial
update supports improvement; the 3e-6 update increases error.

Intervals are conditional on measured association estimates, the fitted models,
and one split and seed. They do not include source effect-estimation error,
refitting, recipe selection, or seed variation. Zero-adapter output agreement
has been checked on a separate bounded sequence set. Native track-output drift
does not measure accuracy against native assay measurements. The complete
native 1-Mb comparison, nested five-seed finalists, and independent external
generalization remain unresolved. These are candidate figures; no adopted
Resource value or figure is changed.
""")


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source",type=Path,required=True)
    parser.add_argument("--check",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    main(parser.parse_args())
