#!/usr/bin/env python3
"""Candidate Fig 6B, 6C, 6D, deep-review fixes (2026-10-17).

Copy of figures/candidates/review-r3-fixes-20261015T112148Z/F6/repair/
panels_20261015T173321Z/f6_render_panels_repair.py (the producer of the installed
6B-6D PDFs), edited to the redraw specs in
Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/PANEL_SPECS.md (6B, 6C, 6D) and
recommended decision 4 of REPORT.md.

What changed relative to the copied producer:
  * Gene calls and classes are read from the adopted Catalog r9d gene table
    (r8 call-state columns, unchanged r8 classes). The r8 reviewed passport
    bundle used by the old producer is read only to assert that the matrix is
    identical.
  * The genetics axis is labelled "eQTL colocalization". The 31 genes with a
    fine-mapped coding variant at PIP >= 0.5 (r9d H87 field coding_cs_state) are
    counted in each 6B cell and shown in a separate 6C column. Classes unchanged.
  * 6B: class rules, units and margins printed; n printed; 0 cell marked as
    zero by construction; tested and untested fills separated (CIE76 dE >= 15);
    footnote for the 7 one-study SuSiE-untestable genes (kept untestable).
  * 6C: SuSiE and ABF in separate columns with study written out; number of
    Tier-1/2 studies with SuSiE PP.H4 > 0.5; RNA drawn as log2FC with 95% CI on
    a shared axis with the 0.5 floor; q printed; routes in the 6D strings.
  * 6D: untested segments in grays, tested classes in one hue; hatched
    sub-segments for the 708 untested genes with one supported call; n printed;
    note on the 20 unresolved genes with a fine-mapped coding variant.

KEY MESSAGE 6B: eQTL colocalization and bulk RNA calls are kept separate; the
three canonical coding MASLD genes sit in "unresolved" because the genetic call
is eQTL colocalization only.
KEY MESSAGE 6C: four worked entries with native estimates; MTARC1 and THRB are
RNA-indeterminate because their change lies inside the 0.5 floor.
KEY MESSAGE 6D: genes routed per class; counts are rule outputs, not priority.

Usage: f6_render_panels_deepfix.py <out_panel_dir> <out_source_table_dir>
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle

FS = 6
plt.rcParams.update({"font.size": FS, "axes.labelsize": FS, "xtick.labelsize": FS, "ytick.labelsize": FS,
                     "pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "Helvetica",
                     "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
                     "hatch.linewidth": 0.4, "text.color": "black", "axes.labelcolor": "black",
                     "xtick.color": "black", "ytick.color": "black"})
OUTLINE = "#757575"      # cell and bar outlines (not the reserved control gray #9E9E9E)
TESTED = "#A9CBEE"       # one tint for every tested class (6B cells, 6D bars)
UNTESTED_6B = "#E6E6E5"  # gray_gradient light end; dE 24.3 / 25.7 / 22.2 vs TESTED (normal/deutan/protan)
UNT_NO_RNA = "#C4C4C4"   # 6D untested, no RNA call
UNT_NO_GEN = "#E6E6E5"   # 6D untested, no genetic call
GE = "≥"            # greater-than-or-equal sign


def mixed_text(ax, xc, y, parts):
    """Centre a line made of (text, fontstyle) parts at data (xc, y)."""
    r = ax.figure.canvas.get_renderer()
    ts = [ax.text(0, y, t, ha="left", va="center", fontstyle=st) for t, st in parts]
    ws = [t.get_window_extent(r).width for t in ts]
    xd = ax.transData.transform((xc, y))[0] - sum(ws) / 2
    for t, w in zip(ts, ws):
        t.set_x(ax.transData.inverted().transform((xd, y))[0])
        xd += w

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
R9D = ROOT / "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/catalog_r9d-20261016T105550Z/build_out/catalog_r9d_gene_table.parquet"
CAT_R8 = ROOT / "figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8"
CATALOG = ROOT / "RNA-seq/results/evidence_passports/releases/gene-catalog-eligibility-20261106"
REC = ROOT / "Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig6/recompute"
RULEBOOK = ROOT / "RNA-seq/results/evidence_passports/candidates/gene-catalog-corrected-r8-20261009T142027Z/passport_rulebook.tsv"
TIERS = ROOT / "GWAS/finemapping/config/gwas_trait_tier.tsv"
OUT = Path(sys.argv[1])
SRC = Path(sys.argv[2])

d = pd.read_parquet(R9D)
assert len(d) == 28233 and d.ensembl_id.is_unique
# The adopted per-gene eligibility and route calls replace only their historical
# equivalents; the r9d assay estimates and coding annotations remain unchanged.
index = pd.read_parquet(CATALOG / "passport_gene_index.parquet").set_index("ensembl_id")
evidence = pd.read_parquet(CATALOG / "passport_evidence_long.parquet")
gen = evidence[evidence.evidence_domain == "genetics"].set_index("ensembl_id")
rna = evidence[evidence.evidence_domain == "transcriptomics"].set_index("ensembl_id")
assert set(index.index) == set(d.ensembl_id) and index.index.is_unique
assert gen.index.is_unique and rna.index.is_unique
assert (d.set_index("ensembl_id").r8_rna_call_state == rna.call_state.reindex(d.ensembl_id)).all()
d["r8_genetics_call_state"] = d.ensembl_id.map(gen.call_state)
d["primary_evidence_class"] = d.ensembl_id.map(index.primary_evidence_class)
routes = {"untested": "EXP_DIRECT_ASSAY_V1", "unresolved": "EXP_UNRESOLVED_V1",
          "established_state_associated": "EXP_STRESS_PERTURB_V1",
          "genetically_anchored": "EXP_ALLELE_AWARE_V1", "concordant": "EXP_MULTI_LAYER_V1"}
d["next_experiment_rule_id"] = d.primary_evidence_class.map(routes)
assert d.primary_evidence_class.value_counts().to_dict() == {
    "untested": 13731, "unresolved": 12958, "established_state_associated": 1118,
    "genetically_anchored": 400, "concordant": 26}
d["coding"] = d.coding_cs_state == "fine_mapped_coding_pip_ge_0.5"
N = len(d)

# cross-check against the r8 reviewed bundle the old producer read
e8 = pd.read_parquet(CAT_R8 / "passport_evidence_long.parquet")
g8 = e8[e8.evidence_domain == "genetics"].set_index("ensembl_id").call_state
r8 = e8[e8.evidence_domain == "transcriptomics"].set_index("ensembl_id").call_state
x = d.set_index("ensembl_id")
assert ((g8.reindex(x.index) == "supported") == (x.r8_genetics_call_state == "supported")).all()
assert (r8.reindex(x.index) == x.r8_rna_call_state).all()

tiers = pd.read_csv(TIERS, sep="\t")
n_t1 = int((tiers.tier == 1).sum()); n_t2 = int((tiers.tier == 2).sum())
assert (n_t1, n_t2) == (15, 20)
units = d[["rna_n_units", "rna_n_control_units", "rna_n_disease_units"]].dropna().drop_duplicates()
assert len(units) == 1
n_units, n_ctl, n_dis = (int(v) for v in units.iloc[0])

states = ["supported", "indeterminate", "untestable"]
G, Rn = d.r8_genetics_call_state, d.r8_rna_call_state
m = pd.crosstab(G, Rn).reindex(index=states, columns=states, fill_value=0)
mc = pd.crosstab(G[d.coding], Rn[d.coding]).reindex(index=states, columns=states, fill_value=0)
assert int(mc.values.sum()) == 31
n_noc = int((d.gen_untestable_reason_r9 == "no_tier12_coloc_record").sum())
n_ins = int((G == "untestable").sum()) - n_noc
assert (n_noc, n_ins) == (9312, 195)
n_below_floor = int(((Rn == "indeterminate") & (d.rna_bh_q < 0.05) & (d.rna_log2fc.abs() <= 0.5)).sum())
assert n_below_floor == 11301

# ---------------- 6B: 3 x 3 call matrix ----------------
label = {("supported", "supported"): "both supported",
         ("supported", "indeterminate"): "genetically\nanchored",
         ("supported", "untestable"): "untested\n(one call supported)",
         ("indeterminate", "supported"): "established\nstate associated",
         ("indeterminate", "indeterminate"): "unresolved",
         ("indeterminate", "untestable"): "untested",
         ("untestable", "supported"): "untested\n(one call supported)",
         ("untestable", "indeterminate"): "untested",
         ("untestable", "untestable"): "untested\nboth untestable"}
EXAMPLES = ["PNPLA3", "TM6SF2", "MTARC1"]
W6B, H6B = 5.1, 3.0
CX0, CW, CH = 1.72, 1.12, 0.60    # cell grid left, cell width, cell height (inches)
ROW_Y = {"supported": 1.73, "indeterminate": 1.13, "untestable": 0.53}   # cell bottoms
row_hdr = {"supported": f"supported\n{int(m.loc['supported'].sum()):,} genes",
           "indeterminate": f"indeterminate\n{int(m.loc['indeterminate'].sum()):,} genes",
           "untestable": f"untestable\n{int(m.loc['untestable'].sum()):,} genes\nno Tier-1/2 COLOC record: {n_noc:,}\nSuSiE untestable in ≥1 study: {n_ins:,}*"}
col_hdr = {"supported": f"supported\n{int(m['supported'].sum()):,} genes",
           "indeterminate": f"indeterminate\n{int(m['indeterminate'].sum()):,} genes\n{n_below_floor:,} with q < 0.05\nbelow the 0.5 floor",
           "untestable": f"untestable\n{int(m['untestable'].sum()):,} genes"}
rows = []
fig = plt.figure(figsize=(W6B, H6B))
ax = fig.add_axes([0, 0, 1, 1]); ax.set(xlim=(0, W6B), ylim=(0, H6B)); ax.axis("off")
for gs in states:
    y = ROW_Y[gs]
    for j, rs in enumerate(states):
        n = int(m.loc[gs, rs]); nc = int(mc.loc[gs, rs])
        tested = gs != "untestable" and rs != "untestable"
        x0 = CX0 + j * CW
        ax.add_patch(Rectangle((x0 + 0.02, y + 0.02), CW - 0.04, CH - 0.04, facecolor=TESTED if tested else UNTESTED_6B,
                               edgecolor=OUTLINE, linewidth=0.4))
        ax.text(x0 + CW / 2, y + 0.47, f"{n:,}", ha="center", va="center")
        ax.text(x0 + CW / 2, y + 0.33, label[(gs, rs)], ha="center", va="center", linespacing=1.0)
        cod_genes = sorted(d.loc[d.coding & (G == gs) & (Rn == rs), "symbol"])
        if nc:
            ex = ", e.g." if (gs, rs) == ("indeterminate", "indeterminate") else ""
            ax.text(x0 + CW / 2, y + 0.17, f"coding PIP {GE} 0.5: {nc}{ex}", ha="center", va="center")
        if gs == "indeterminate" and rs == "indeterminate":
            assert all(s in cod_genes for s in EXAMPLES)
            ax.text(x0 + CW / 2, y + 0.085, ", ".join(EXAMPLES), ha="center", va="center", fontstyle="italic")
        rows.append(dict(record="cell", genetics_call=gs, rna_call=rs, n_genes=n, n_coding_pip_ge_0_5=nc,
                         coding_genes=";".join(cod_genes), display_label=label[(gs, rs)].replace("\n", " "),
                         catalog_class=";".join(sorted(d.loc[(G == gs) & (Rn == rs), "primary_evidence_class"].astype(str).unique())),
                         tested_both=tested, fill=TESTED if tested else UNTESTED_6B))
    ax.text(CX0 - 0.05, y + CH / 2, row_hdr[gs], ha="right", va="center", linespacing=1.0)
for j, rs in enumerate(states):
    ax.text(CX0 + j * CW + CW / 2, 2.34, col_hdr[rs], ha="center", va="bottom", linespacing=1.0)
ax.text(CX0 + 1.5 * CW, 2.87, f"Bulk liver RNA, disease vs control: {n_units} genetic-identity clusters\n"
        f"({n_ctl} control, {n_dis} disease); supported = BH q < 0.05 and |log2FC| > 0.5",
        ha="center", va="center", linespacing=1.0)
ax.text(0.17, 1.43, f"Liver-trait eQTL colocalization (SuSiE-COLOC)\n{n_t1 + n_t2} Tier-1/2 GWAS: {n_t1} MASLD, {n_t2} liver enzyme\n"
        "supported = PP.H4 > 0.5", ha="center", va="center", rotation=90, linespacing=1.0)
ax.text(0.05, 0.40, f"n = {N:,} genes. Coding PIP {GE} 0.5: fine-mapped coding variant at SuSiE PIP {GE} 0.5 "
        f"(Catalog r9d, {int(d.coding.sum())} genes),", ha="left", va="center")
ax.text(0.05, 0.30, "not part of the genetic call.", ha="left", va="center")
ax.text(0.05, 0.17, "* An untestable study makes genetics untestable unless another study provides positive support.",
        ha="left", va="center")
fig.savefig(OUT / "fig6b_coverage_observability.pdf")
