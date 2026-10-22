#!/usr/bin/env python3
"""Candidate Fig 6B, 6C, 6D for review round 3 (workstream F6).

Inputs: Catalog r8 reviewed copy (unchanged calls), the adopted 774-unit bulk
coefficient table and COLOC r2 gene-level table, via the F6 re-derivation
outputs. No counts are changed except where stated in the source tables.
"""
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle

plt.rcParams.update({"font.size": 6, "axes.labelsize": 6, "xtick.labelsize": 6, "ytick.labelsize": 6,
                     "pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "Helvetica",
                     "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5})
GRAY = "#9E9E9E"
LIGHT = "#E8EEF5"
LGRAY = "#EEEEEE"
ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
CAT = ROOT / "figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8"
RED = Path(sys.argv[1])   # rederive dir
OUT = Path(sys.argv[2])   # panel output dir (must exist and be empty of these files)

e = pd.read_parquet(CAT / "passport_evidence_long.parquet")
g = pd.read_parquet(CAT / "passport_gene_index.parquet").set_index("ensembl_id")
gen = e[e.evidence_domain == "genetics"].set_index("ensembl_id")
rna = e[e.evidence_domain == "transcriptomics"].set_index("ensembl_id")
side = pd.read_csv(RED / "R15_03_gene_effect_sidecar.tsv.gz", sep="\t").set_index("ensembl_id")

# ---------------- 6B: 3 x 3 call matrix, no tested-negative row/column ----------------
states = ["supported", "indeterminate", "untestable"]
m = pd.crosstab(gen.call_state, rna.call_state).reindex(index=states, columns=states, fill_value=0)
label = {("supported", "supported"): "both supported",
         ("supported", "indeterminate"): "genetically\nanchored",
         ("supported", "untestable"): "untested\n(one call supported)",
         ("indeterminate", "supported"): "established\nstate associated",
         ("indeterminate", "indeterminate"): "unresolved",
         ("indeterminate", "untestable"): "untested",
         ("untestable", "supported"): "untested\n(one call supported)",
         ("untestable", "indeterminate"): "untested",
         ("untestable", "untestable"): "untested"}
rows = []
for gs in states:
    for rs in states:
        rows.append(dict(genetics_call=gs, rna_call=rs, n_genes=int(m.loc[gs, rs]), display_label=label[(gs, rs)].replace("\n", " "),
                         catalog_class=g.loc[gen.index[(gen.call_state == gs)].intersection(rna.index[rna.call_state == rs]), "primary_evidence_class"].unique().tolist()))
src6b = pd.DataFrame(rows)
src6b.to_csv(OUT / "fig6b_call_matrix_3x3_source.tsv", sep="\t", index=False)
col_hdr = {"supported": "supported", "indeterminate": "indeterminate", "untestable": "untestable\n(not detected)"}
row_hdr = {"supported": "supported", "indeterminate": "indeterminate", "untestable": "untestable\n(no COLOC record;\n7 SuSiE untestable)"}
fig, ax = plt.subplots(figsize=(3.8, 2.6))
ax.set(xlim=(-1.15, 3), ylim=(0, 3.45)); ax.axis("off")
fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
for i, gs in enumerate(states):
    y = 2 - i
    for j, rs in enumerate(states):
        n = int(m.loc[gs, rs])
        inferential = gs != "untestable" and rs != "untestable"
        ax.add_patch(Rectangle((j + 0.03, y + 0.03), 0.94, 0.94, facecolor=LIGHT if inferential else LGRAY, edgecolor=GRAY, linewidth=0.4))
        ax.text(j + 0.5, y + 0.62, f"{n:,}", ha="center", va="center")
        ax.text(j + 0.5, y + 0.36, label[(gs, rs)], ha="center", va="center", linespacing=1.0)
    ax.text(-0.05, y + 0.5, row_hdr[gs], ha="right", va="center", linespacing=1.0)
for j, rs in enumerate(states):
    ax.text(j + 0.5, 3.1, col_hdr[rs], ha="center", va="center", linespacing=1.0)
ax.text(1.5, 3.38, "Established disease-state transcriptome (assay call)", ha="center", va="center")
ax.text(-1.1, 1.5, "Liver-trait regulatory genetics (assay call)", ha="center", va="center", rotation=90)
fig.savefig(OUT / "fig6b_coverage_observability_candidate.pdf")
plt.close(fig)

# ---------------- 6C: worked entries with native estimates ----------------
genes = ["GLP1R", "HKDC1", "MTARC1", "THRB"]
s = side[side.symbol.isin(genes)].copy()


def gtxt(r):
    if pd.notna(r.coloc_best_susie_pp4):
        return f"{r.coloc_best_susie_pp4:.3f} SuSiE\n{r.coloc_best_susie_study}"
    if pd.notna(r.coloc_best_abf_pp4):
        return f"no SuSiE pair;\nABF {r.coloc_best_abf_pp4:.3f}"
    return "no COLOC record"


def qfmt(q):
    return f"{q:.1e}".replace("e-0", "e-") if q < 0.001 else f"{q:.3f}"


def rtxt(r):
    if pd.isna(r.rna_log2fc):
        return "not detected"
    return f"{r.rna_log2fc:+.2f} (q {qfmt(r.rna_bh_q)})"


cls = {"genetically_anchored": "Genetically anchored", "concordant": "Both supported", "unresolved": "Unresolved", "untested": "Untested"}
s["genetics_estimate_text"] = s.apply(gtxt, axis=1).str.replace("\n", " ")
s["rna_estimate_text"] = s.apply(rtxt, axis=1)
s["class_display"] = s.primary_evidence_class.map(cls)
s.reset_index()[["ensembl_id", "symbol", "genetics_call", "coloc_best_susie_pp4", "coloc_best_susie_study", "coloc_best_abf_pp4",
                 "coloc_best_abf_study", "genetics_estimate_text", "rna_call", "rna_log2fc", "rna_se", "rna_bh_q", "rna_estimate_text",
                 "primary_evidence_class", "class_display"]].sort_values("symbol").to_csv(OUT / "fig6c_worked_entries_candidate_source.tsv", sep="\t", index=False)
route = {"GLP1R": "Acquire missing decisive assay\nMeasure missing assay",
         "HKDC1": "Allele edit + gene perturbation\nCoding-state perturbation",
         "MTARC1": "Measure decisive missing layer\nCoding-state perturbation",
         "THRB": "Allele-aware editing\nCoding-state perturbation"}
fig, ax = plt.subplots(figsize=(4.7, 1.9))
ax.set(xlim=(0, 1), ylim=(0, 1)); ax.axis("off")
fig.subplots_adjust(left=0.01, right=0.99, top=0.98, bottom=0.02)
x = [0.0, 0.09, 0.335, 0.555, 0.735, 0.99]
hdr = ["Gene", "Genetics call\nbest PP.H4, method, study", "RNA call\nlog2FC (BH q), 774 units", "Gene evidence class", "Gene route /\ncoding-object route", ""]
for xx, h in zip(x[:5], hdr[:5]):
    ax.text(xx, 0.92, h, ha="left", va="center", linespacing=1.0)
for i, sym in enumerate(sorted(genes)):
    r = s[s.symbol == sym].iloc[0]
    y = 0.72 - i * 0.19
    ax.text(x[0], y, sym, ha="left", va="center", fontstyle="italic")
    ax.text(x[1], y, r.genetics_call + "\n" + gtxt(r), ha="left", va="center", linespacing=1.0)
    ax.text(x[2], y, r.rna_call + "\n" + rtxt(r), ha="left", va="center", linespacing=1.0)
    ax.text(x[3], y, cls[r.primary_evidence_class], ha="left", va="center")
    ax.text(x[4], y, route[sym], ha="left", va="center", linespacing=1.0)
fig.savefig(OUT / "fig6c_worked_entries_candidate.pdf")
plt.close(fig)

# ---------------- 6D: routed counts, untested split by missing domain ----------------
untested = g.index[g.primary_evidence_class == "untested"]
miss_rna = int((rna.loc[untested, "call_state"] == "untestable").sum())
miss_gen = int((gen.loc[untested, "call_state"] == "untestable").sum())
both_missing = int(((rna.loc[untested, "call_state"] == "untestable") & (gen.loc[untested, "call_state"] == "untestable")).sum())
cc = g.primary_evidence_class.value_counts()
bars = [("untested", "acquire the missing decisive assay", [("no RNA call", miss_rna - both_missing), ("no genetic call", miss_gen - both_missing)]),
        ("unresolved", "measure the single most decisive missing layer", [("", int(cc["unresolved"]))]),
        ("established state associated", "gene knockdown and rescue", [("", int(cc["established_state_associated"]))]),
        ("genetically anchored", "allele-aware regulatory editing", [("", int(cc["genetically_anchored"]))]),
        ("both supported", "allele-aware edit plus gene perturbation", [("", int(cc["concordant"]))])]
src6d = []
for name, rule, segs in bars:
    for seg, n in segs:
        src6d.append(dict(primary_evidence_class_display=name, route=rule, segment=seg or name, n_genes=n))
src6d = pd.DataFrame(src6d)
src6d.loc[len(src6d)] = dict(primary_evidence_class_display="check", route="untested genes with both domains untestable", segment="both missing", n_genes=both_missing)
src6d.to_csv(OUT / "fig6d_next_experiment_candidate_source.tsv", sep="\t", index=False)
fig, ax = plt.subplots(figsize=(4.6, 2.2))
fig.subplots_adjust(left=0.42, right=0.95, top=0.97, bottom=0.17)
ypos = np.arange(len(bars))[::-1]
for yy, (name, rule, segs) in zip(ypos, bars):
    left = 0
    for k, (seg, n) in enumerate(segs):
        ax.barh(yy, n, left=left, height=0.6, color="#C9D4E0" if k == 0 else "#7F95AD", edgecolor="black", linewidth=0.3)
        if seg:
            ax.text(left + n / 2, yy, f"{seg}\n{n:,}", ha="center", va="center", linespacing=1.0)
        left += n
    ax.text(left + 150, yy, f"{left:,}", ha="left", va="center")
ax.set_yticks(ypos)
ax.set_yticklabels([f"{n}\n{r}" for n, r, _ in bars], linespacing=1.0)
ax.set_xlabel("Genes routed (count, not priority)")
ax.set_xlim(0, 16000)
ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, p: f"{v:,.0f}"))
ax.set_xticks([0, 5000, 10000, 15000])
ax.spines[["top", "right"]].set_visible(False)
ax.tick_params(axis="y", length=0)
fig.savefig(OUT / "fig6d_next_experiment_candidate.pdf")
plt.close(fig)
print("6B", m.to_dict(), "6D untested split", miss_rna, miss_gen, both_missing)
