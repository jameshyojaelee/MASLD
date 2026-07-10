#!/usr/bin/env python
"""plot_offtarget_report.py — guide QC report for the Cas13 vM38 library.

Only the checks that are actually meaningful:
  off-target — the raw genome/transcriptome off-target DISTRIBUTIONS (genome-level is
    informational-only; the transcriptome screen is the hard filter). Per-category
    "genome-clean rate" breakdowns are deliberately omitted: genome perfect-matches
    live overwhelmingly in introns/intergenic/repeats that an RNA-targeting enzyme
    never sees, so "genome-clean %" is not a guide-quality axis.
  on-target — isoform coverage.
  efficacy  — TIGER + Cas13Design score distributions.

Data: data/guides/cache/offtarget_merged.parquet (built by offtarget_tx_vs_genome.py).
House style: all text black; marks may be colored; bars carry counts/percentages; no lollipops.
Output: figures/guide_qc/offtarget_check_report.pdf
"""
import polars as pl, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Cas13_Library_Design"
PARQ = f"{ROOT}/data/guides/cache/offtarget_merged.parquet"
OUT  = f"{ROOT}/figures/guide_qc/offtarget_check_report.pdf"

CLEAN="#4C9F70"; MM0="#C0143C"; MM1="#E88A2A"; MM2="#F2C744"; INDEL="#7F7F7F"; PC_C="#3B6FB6"

d = pl.read_parquet(PARQ).with_columns([
    pl.col("g0").fill_null(0), pl.col("g1").fill_null(0), pl.col("g2").fill_null(0)])
d = d.with_columns(((pl.col("g0")==0)&(pl.col("g1")==0)&(pl.col("g2")==0)).alias("clean"))
pdf = d.to_pandas()
N = len(pdf)

def barlab(ax, bars, vals, fs=8):
    for b,v in zip(bars, vals):
        ax.text(b.get_x()+b.get_width()/2, b.get_height(),
                f"{int(v)}\n{100*v/N:.1f}%", ha="center", va="bottom", fontsize=fs)

fig, axs = plt.subplots(2, 4, figsize=(19, 9.6))
fig.subplots_adjust(hspace=0.42, wspace=0.28, top=0.86, bottom=0.09, left=0.05, right=0.99)

# ============================ row 1: OFF-TARGET ============================
# P1: genome-level off-target tiers
ax=axs[0,0]
cats=["0:0:0\n(no OT)","≥1 @0mm","≥1 @1mm","≥1 @2mm","any\nindel"]
indel=int((pdf.genome_any_indel.astype(str).str.lower()=="true").sum())
vals=[int(pdf.clean.sum()), int((pdf.g0>0).sum()), int((pdf.g1>0).sum()), int((pdf.g2>0).sum()), indel]
b=ax.bar(cats, vals, color=[CLEAN,MM0,MM1,MM2,INDEL]); barlab(ax,b,vals,fs=7)
ax.set_title("Genome-level off-target tiers", fontsize=11, fontweight="bold")
ax.set_ylabel("guides"); ax.set_ylim(0, N*1.15)

# P2: genome vs transcriptome (why genome is inflated)
ax=axs[0,1]
g_any=int((pdf.g0>0).sum())
tx_any=int(((pdf.tx_offtarget_pseudogene_only.astype(str).str.lower()=="true")
            | (pdf.tx_gene_family_hit_transcripts.fillna("")!="")).sum())
b=ax.bar(["Genome\n(informational)","Transcriptome\n(hard filter)"], [g_any, tx_any],
         color=[MM0, CLEAN]); barlab(ax,b,[g_any,tx_any])
ax.set_title("≥1 perfect-match off-target:\ngenome vs transcriptome", fontsize=11, fontweight="bold")
ax.set_ylabel("guides"); ax.set_ylim(0, N*1.15)

# P3: transcriptome off-target status
ax=axs[0,2]
tx_pseudo=int((pdf.tx_offtarget_pseudogene_only.astype(str).str.lower()=="true").sum())
tx_gf=int(((pdf.tx_gene_family_hit_transcripts.fillna("")!="")|(pdf.tx_gene_family_hit_genes.fillna("")!="")).sum())
tx_clean=N - tx_pseudo - tx_gf
cats=["fully\nclean","gene-family\n(tolerated)","pseudogene\n(tolerated)","disquali-\nfying"]
vals=[tx_clean, tx_gf, tx_pseudo, 0]
b=ax.bar(cats, vals, color=[CLEAN,MM1,MM2,MM0]); barlab(ax,b,vals,fs=7)
ax.set_title("Transcriptome off-target status", fontsize=11, fontweight="bold")
ax.set_ylabel("guides"); ax.set_ylim(0, N*1.15)

# P4: genome 0mm off-target count per guide (raw distribution)
ax=axs[0,3]
g0=pdf.g0.clip(upper=6)
cnt,_=np.histogram(g0,bins=np.arange(0,8)-0.5)
xt=[str(i) for i in range(6)]+["6+"]
b=ax.bar(range(7), cnt, color=[CLEAN]+[MM0]*6)
for bar,v in zip(b,cnt):
    if v>0: ax.text(bar.get_x()+bar.get_width()/2,v,f"{int(v)}",ha="center",va="bottom",fontsize=7)
ax.set_xticks(range(7)); ax.set_xticklabels(xt)
ax.set_title("Genome 0mm off-target count / guide", fontsize=11, fontweight="bold")
ax.set_xlabel("# 0mm genome off-targets"); ax.set_ylabel("guides")

# ============================ row 2: ON-TARGET + EFFICACY ============================
# P5: on-target isoform coverage
ax=axs[1,0]
nt=pdf.n_target.clip(upper=9)
cnt,_=np.histogram(nt,bins=np.arange(1,11)-0.5)
xt=[str(i) for i in range(1,9)]+["9+"]
ax.bar(range(1,10), cnt, color=PC_C)
ax.set_xticks(range(1,10)); ax.set_xticklabels(xt)
ax.set_title("On-target isoform coverage (n_target)", fontsize=11, fontweight="bold")
ax.set_xlabel("# target isoforms hit"); ax.set_ylabel("guides")

# P6: TIGER score
ax=axs[1,1]
tg=pdf.tiger_score.dropna()
ax.hist(tg, bins=np.linspace(0,1,41), color=PC_C, edgecolor="white")
ax.axvline(0.75, color=MM0, ls="--", lw=1.2)
ax.text(0.74, ax.get_ylim()[1]*0.92, "≥0.75 gate ", color="black", fontsize=7.5, ha="right")
ax.set_title(f"TIGER score (median {tg.median():.3f})", fontsize=11, fontweight="bold")
ax.set_xlabel("TIGER score"); ax.set_ylabel("guides"); ax.set_xlim(0,1)

# P7: Cas13Design score
ax=axs[1,2]
cd=pdf.cas13_score.dropna()
ax.hist(cd, bins=np.linspace(0,1,41), color=CLEAN, edgecolor="white")
ax.axvline(0.75, color=MM0, ls="--", lw=1.2)
ax.text(0.74, ax.get_ylim()[1]*0.92, "≥0.75 gate ", color="black", fontsize=7.5, ha="right")
ax.set_title(f"Cas13Design score (median {cd.median():.3f})", fontsize=11, fontweight="bold")
ax.set_xlabel("Cas13Design score"); ax.set_ylabel("guides"); ax.set_xlim(0,1)

# P8: TIGER vs Cas13Design (joint; OR-gate)
ax=axs[1,3]
j=pdf.dropna(subset=["tiger_score","cas13_score"])
hb=ax.hexbin(j.tiger_score, j.cas13_score, gridsize=40, cmap="viridis", mincnt=1, bins="log")
ax.axvline(0.75, color=MM0, ls="--", lw=1.0); ax.axhline(0.75, color=MM0, ls="--", lw=1.0)
ax.set_title("TIGER vs Cas13Design (gate: either ≥0.75)", fontsize=11, fontweight="bold")
ax.set_xlabel("TIGER score"); ax.set_ylabel("Cas13Design score"); ax.set_xlim(0,1); ax.set_ylim(0,1)
cb=fig.colorbar(hb, ax=ax, fraction=0.046, pad=0.04); cb.set_label("guides (log)", fontsize=7)

fig.suptitle(f"Cas13 vM38 library — guide QC report  ({N:,} guides, "
             f"{pdf.gene_symbol_mouse.nunique():,} genes)", fontsize=14, fontweight="bold", y=0.975)
fig.savefig(OUT, bbox_inches="tight"); plt.close(fig)
print("saved ->", OUT)
