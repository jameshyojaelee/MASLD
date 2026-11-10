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
plt.close(fig)
for gs in states:
    rows.append(dict(record="row_margin", genetics_call=gs, rna_call="All", n_genes=int(m.loc[gs].sum()), n_coding_pip_ge_0_5=int(mc.loc[gs].sum())))
for rs in states:
    rows.append(dict(record="col_margin", genetics_call="All", rna_call=rs, n_genes=int(m[rs].sum()), n_coding_pip_ge_0_5=int(mc[rs].sum())))
rows += [dict(record="total", genetics_call="All", rna_call="All", n_genes=N, n_coding_pip_ge_0_5=int(d.coding.sum())),
         dict(record="annotation", genetics_call="untestable", rna_call="All", display_label="no Tier-1/2 COLOC record", n_genes=n_noc),
         dict(record="annotation", genetics_call="untestable", rna_call="All", display_label="SuSiE untestable in at least one Tier-1/2 study", n_genes=n_ins),
         dict(record="annotation", genetics_call="All", rna_call="indeterminate", display_label="BH q < 0.05 and |log2FC| <= 0.5", n_genes=n_below_floor),
         dict(record="annotation", genetics_call="All", rna_call="All", display_label=f"bulk units {n_units} = {n_ctl} control + {n_dis} disease genetic-identity clusters", n_genes=n_units),
         dict(record="annotation", genetics_call="All", rna_call="All", display_label=f"Tier-1/2 GWAS {n_t1 + n_t2} = {n_t1} direct MASLD + {n_t2} liver enzyme", n_genes=n_t1 + n_t2)]
pd.DataFrame(rows).to_csv(SRC / "6B.tsv", sep="\t", index=False)

# ---------------- route strings (shared by 6C and 6D) ----------------
rb = pd.read_csv(RULEBOOK, sep="\t").set_index("rule_id").perturbation
GENE_ROUTE = {"untested": "acquire the missing decisive assay", "unresolved": "measure the single most decisive missing layer",
              "established_state_associated": "gene knockdown and rescue", "genetically_anchored": "allele-aware regulatory editing",
              "concordant": "allele-aware edit plus gene perturbation"}
CLASS_RULE = {"untested": "EXP_DIRECT_ASSAY_V1", "unresolved": "EXP_UNRESOLVED_V1", "established_state_associated": "EXP_STRESS_PERTURB_V1",
              "genetically_anchored": "EXP_ALLELE_AWARE_V1", "concordant": "EXP_MULTI_LAYER_V1"}
for c, rid in CLASS_RULE.items():
    assert rb[rid].startswith(GENE_ROUTE[c]), (rid, rb[rid])
    assert set(d.loc[d.primary_evidence_class == c, "next_experiment_rule_id"]) == {rid}
OBJ_ROUTE = {"protein_coding_gene": "gene perturbation + rescue", "untestable_nomination": "acquire the missing decisive measurement"}

# ---------------- 6C: worked entries ----------------
genes = ["GLP1R", "HKDC1", "MTARC1", "THRB"]
ps = pd.read_csv(REC / "fig6c_worked_entries_per_study.tsv", sep="\t")
ps = ps[ps.gwas_name.isin(tiers.loc[tiers.tier.isin([1, 2]), "study_name"])]
obj = pd.read_csv(REC / "fig6c_object_routes.tsv", sep="\t").set_index("symbol")
STUDY = {"UKBB_AST": "UK Biobank AST", "UKBB_GGT": "UK Biobank GGT", "UKBB_ALT": "UK Biobank ALT", "MVP_NAFLD_EAS": "MVP NAFLD, EAS",
         "FinnGen_NASH": "FinnGen NASH", "2022_36402844_PDFF_EUR": "PDFF GWAS"}


def pfmt(v):
    return "> 0.999" if v > 0.9995 else f"{v:.3f}"


def qfmt(q):
    return f"{q:.1e}".replace("e-0", "e-") if q < 0.001 else f"{q:.3f}"


CLS = {"genetically_anchored": "Genetically\nanchored", "concordant": "Both supported*", "unresolved": "Unresolved", "untested": "Untested"}
rows6c = []
for sym in genes:
    r = d[d.symbol == sym].iloc[0]
    p = ps[ps.symbol == sym]
    n_su = int(p["PP.H4.susie"].notna().sum()); n_su_pos = int((p["PP.H4.susie"] > 0.5).sum())
    if pd.notna(r.gen_coloc_best_susie_pp4):
        assert np.isclose(p["PP.H4.susie"].max(), r.gen_coloc_best_susie_pp4)
        susie_txt = f"{pfmt(r.gen_coloc_best_susie_pp4)}, {STUDY[r.gen_coloc_best_susie_study]}\n> 0.5 in {n_su_pos} of {n_su} studies"
    else:
        assert n_su == 0
        susie_txt = "no Tier-1/2 SuSiE result"
    assert np.isclose(p["PP.H4.abf"].max(), r.gen_coloc_best_abf_pp4)
    abf_txt = f"{pfmt(r.gen_coloc_best_abf_pp4)}\n{STUDY[r.gen_coloc_best_abf_study]}"
    if r.coding:
        cod_txt = f"{r.coding_cs_top_consequence.replace('_variant', '')}, PIP {r.coding_cs_max_pip:.2f}\n{STUDY[r.coding_cs_top_study]}\nPIP {GE} 0.5 in {int(r.coding_cs_n_studies)} studies"
    else:
        assert r.coding_cs_state == "no_coding_cs_member"
        cod_txt = "none"
    if pd.notna(r.rna_log2fc):
        lo, hi = r.rna_log2fc - 1.96 * r.rna_se, r.rna_log2fc + 1.96 * r.rna_se
        q_txt = f"{r.rna_log2fc:+.2f}\nq {qfmt(r.rna_bh_q)}".replace("+-", "-")
    else:
        lo = hi = np.nan
        q_txt = ""
    o = obj.loc[sym, "candidate_object"]
    rows6c.append(dict(symbol=sym, ensembl_id=r.ensembl_id, genetics_call_eqtl_coloc=r.r8_genetics_call_state,
                       susie_best_pp4=r.gen_coloc_best_susie_pp4, susie_best_study=r.gen_coloc_best_susie_study,
                       n_tier12_studies_with_susie=n_su, n_tier12_studies_susie_gt_0_5=n_su_pos,
                       abf_best_pp4=r.gen_coloc_best_abf_pp4, abf_best_study=r.gen_coloc_best_abf_study,
                       coding_cs_state=r.coding_cs_state, coding_max_pip=r.coding_cs_max_pip, coding_top_study=r.coding_cs_top_study,
                       coding_consequence=r.coding_cs_top_consequence, coding_n_studies_pip_ge_0_5=r.coding_cs_n_studies,
                       rna_call=r.r8_rna_call_state, rna_log2fc=r.rna_log2fc, rna_se=r.rna_se, rna_ci95_lo=lo, rna_ci95_hi=hi,
                       rna_bh_q=r.rna_bh_q, primary_evidence_class=r.primary_evidence_class, class_display=CLS[r.primary_evidence_class].replace("\n", " "),
                       gene_route_rule=r.next_experiment_rule_id, gene_route=GENE_ROUTE[r.primary_evidence_class],
                       molecular_object=o, object_route_rule=obj.loc[sym, "next_experiment_rule_id"], object_route=OBJ_ROUTE[o],
                       susie_text=susie_txt.replace("\n", " | "), abf_text=abf_txt.replace("\n", " | "), coding_text=cod_txt.replace("\n", " | "),
                       rna_text=q_txt.replace("\n", " | ") if q_txt else "not detected"))
c6 = pd.DataFrame(rows6c)
c6.to_csv(SRC / "6C.tsv", sep="\t", index=False)

W6C, H6C = 5.5, 2.7
XG, XS, XA, XC, XQ, XK = 0.04, 0.46, 1.48, 2.30, 4.36, 4.82
PX0, PX1, VMIN, VMAX = 3.28, 4.28, -0.6, 1.6
fx = lambda v: PX0 + (v - VMIN) / (VMAX - VMIN) * (PX1 - PX0)
ROWTOP, RH, AXY = 2.33, 0.44, 0.55
fig = plt.figure(figsize=(W6C, H6C))
ax = fig.add_axes([0, 0, 1, 1]); ax.set(xlim=(0, W6C), ylim=(0, H6C)); ax.axis("off")
ax.text(XS, 2.62, f"eQTL colocalization, {n_t1 + n_t2} Tier-1/2 GWAS", ha="left", va="center")
ax.plot([XS, XC - 0.08], [2.56, 2.56], color="black", lw=0.4)
ax.text(XG, 2.45, "Gene", ha="left", va="center")
ax.text(XS, 2.45, "SuSiE PP.H4, study", ha="left", va="center")
ax.text(XA, 2.45, "ABF PP.H4, study", ha="left", va="center")
ax.text(XC, 2.50, "Fine-mapped coding\nvariant", ha="left", va="center", linespacing=1.0)
ax.text(PX0, 2.53, f"Bulk liver RNA\n{n_units} genetic-identity clusters\n({n_ctl} control, {n_dis} disease)", ha="left", va="center", linespacing=1.0)
ax.text(XK, 2.45, "Class", ha="left", va="center")
yc = {}
for i, sym in enumerate(genes):
    r = c6[c6.symbol == sym].iloc[0]
    top = ROWTOP - i * RH
    y = top - 0.12; yc[sym] = y
    ax.text(XG, y, sym, ha="left", va="center", fontstyle="italic")
    ax.text(XS, y, r.susie_text.replace(" | ", "\n"), ha="left", va="center", linespacing=1.0)
    ax.text(XA, y, r.abf_text.replace(" | ", "\n"), ha="left", va="center", linespacing=1.0)
    ax.text(XC, y, r.coding_text.replace(" | ", "\n"), ha="left", va="center", linespacing=1.0)
    ax.text(XK, y, CLS[r.primary_evidence_class], ha="left", va="center", linespacing=1.0)
    ax.text(XS, top - 0.36, f"Route: gene level, {r.gene_route}; molecular object, {r.object_route}", ha="left", va="center")
    if np.isfinite(r.rna_log2fc):
        ax.plot([fx(r.rna_ci95_lo), fx(r.rna_ci95_hi)], [y, y], color="black", lw=0.6, solid_capstyle="butt")
        ax.plot([fx(r.rna_log2fc)], [y], marker="o", ms=2.6, mfc="black", mec="black", mew=0, ls="none")
        ax.text(XQ, y, r.rna_text.replace(" | ", "\n"), ha="left", va="center", linespacing=1.0)
    else:
        ax.text((PX0 + PX1) / 2, y, "not detected", ha="center", va="center")
    if i:
        ax.plot([XG, W6C - 0.04], [top + 0.005, top + 0.005], color="#D0D0D0", lw=0.3)
for v, ls, col in [(-0.5, (0, (2, 1.5)), "black"), (0.5, (0, (2, 1.5)), "black"), (0.0, "-", "#BDBDBD")]:
    for sym in genes[1:]:   # reference lines only inside the data line of rows with an RNA estimate
        ax.plot([fx(v), fx(v)], [yc[sym] - 0.14, yc[sym] + 0.14], color=col, lw=0.4, ls=ls, zorder=0)
    ax.plot([fx(v), fx(v)], [AXY, AXY + 0.06], color=col, lw=0.4, ls=ls, zorder=0)
ax.plot([PX0, PX1], [AXY, AXY], color="black", lw=0.5)
TICKS = [-0.5, 0.0, 0.5, 1.0, 1.5]
for t in TICKS:
    ax.plot([fx(t), fx(t)], [AXY, AXY - 0.035], color="black", lw=0.5)
    ax.text(fx(t), AXY - 0.06, f"{t:g}".replace("-", "−"), ha="center", va="top")
ax.text((PX0 + PX1) / 2, AXY - 0.18, "log2FC, disease vs control (95% CI)", ha="center", va="center")
ax.text(XG, 0.24, "RNA supported = BH q < 0.05 and |log2FC| > 0.5 (dashed lines). Genetic call = SuSiE PP.H4 > 0.5; ABF for reference only.",
        ha="left", va="center")
ax.text(XG, 0.14, "Coding variant: Catalog r9d field, not part of the genetic call.", ha="left", va="center")
ax.text(XG, 0.04, "* Requires a supported call in each domain; direction not compared.", ha="left", va="center")
fig.savefig(OUT / "fig6c_worked_entries.pdf")
plt.close(fig)
pd.DataFrame({"tick_log2fc": TICKS, "x_in": [fx(t) for t in TICKS]}).to_csv(SRC / "6C_axis.tsv", sep="\t", index=False)

# ---------------- 6D: routed counts ----------------
cc = d.primary_evidence_class.value_counts()
unt = d.primary_evidence_class == "untested"
no_rna = unt & (Rn == "untestable")
no_gen = unt & (G == "untestable") & (Rn != "untestable")
assert not (no_rna & no_gen).any()
both_untestable = int(((Rn == "untestable") & (G == "untestable")).sum())
n_no_rna, n_no_gen = int(no_rna.sum()), int(no_gen.sum())
one_gen = int((no_rna & (G == "supported")).sum()); one_rna = int((no_gen & (Rn == "supported")).sum())
assert (n_no_rna, n_no_gen, one_gen, one_rna, both_untestable) == (4232, 9499, 10, 709, 8)
n_unres_coding = int((d.coding & (d.primary_evidence_class == "unresolved")).sum())
assert n_unres_coding == 20
order = ["untested", "unresolved", "established_state_associated", "genetically_anchored", "concordant"]
NAME = {"untested": "untested", "unresolved": "unresolved", "established_state_associated": "established state associated",
        "genetically_anchored": "genetically anchored", "concordant": "both supported"}
segs = {"untested": [("RNA untestable", n_no_rna, UNT_NO_RNA, one_gen, "genetics supported, RNA not detected"),
                     ("genetics untestable only", n_no_gen, UNT_NO_GEN, one_rna, "RNA supported, genetics untestable")]}
for c in order[1:]:
    segs[c] = [("", int(cc[c]), TESTED, 0, "")]
src6d = []
fig = plt.figure(figsize=(5.2, 2.5))
ax = fig.add_axes([0.405, 0.30, 0.565, 0.66])
ax.set_xlim(0, 16000)
ax.set_ylim(-0.5, 4.5)
ypos = np.arange(len(order))[::-1]
for yy, c in zip(ypos, order):
    left = 0
    for seg, n, col, nh, hdesc in segs[c]:
        ax.barh(yy, n, left=left, height=0.7, color=col, edgecolor=OUTLINE, linewidth=0.3)
        src6d.append(dict(primary_evidence_class=c, class_display=NAME[c], gene_route_rule=CLASS_RULE[c], gene_route=GENE_ROUTE[c],
                          segment=seg or NAME[c], n_genes=n, bar_left=left, bar_right=left + n, fill=col, hatched=False))
        if nh:
            ax.barh(yy, nh, left=left, height=0.7, color="none", edgecolor="black", linewidth=0.0, hatch="//////")
            src6d.append(dict(primary_evidence_class=c, class_display=NAME[c], gene_route_rule=CLASS_RULE[c], gene_route=GENE_ROUTE[c],
                              segment=f"one supported call: {hdesc}", n_genes=nh, bar_left=left, bar_right=left + nh, fill="none", hatched=True))
        if seg:
            ax.text(left + n / 2 + (350 if nh > 100 else 0), yy, f"{seg}\n{n:,}", ha="center", va="center", linespacing=1.0)
        left += n
    ax.text(left + 150, yy, f"{left:,}", ha="left", va="center")
    if c == "unresolved":
        ax.text(left / 2, yy + 0.13, f"includes {n_unres_coding} genes with a fine-mapped coding variant:", ha="center", va="center")
        mixed_text(ax, left / 2, yy - 0.15, [("e.g.", "normal"), (" " + ", ".join(EXAMPLES), "italic")])
ax.set_yticks(ypos)
ax.set_yticklabels([f"{NAME[c]}\n{GENE_ROUTE[c]}" for c in order], linespacing=1.0)
ax.set_xlabel(f"Genes routed (count, not priority)\nn = {N:,} genes; genetic call = liver-trait eQTL colocalization", linespacing=1.1)
ax.set_xlim(0, 16000)
ax.set_ylim(-0.5, len(order) - 0.5)
ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, p: f"{v:,.0f}"))
ax.set_xticks([0, 5000, 10000, 15000])
ax.spines[["top", "right"]].set_visible(False)
ax.tick_params(axis="y", length=0)
lax = fig.add_axes([0.405, 0.04, 0.565, 0.07]); lax.axis("off"); lax.set(xlim=(0, 1), ylim=(0, 1))
lax.add_patch(Rectangle((0.0, 0.2), 0.035, 0.6, facecolor="none", edgecolor="black", linewidth=0.0, hatch="//////"))
lax.add_patch(Rectangle((0.0, 0.2), 0.035, 0.6, facecolor="none", edgecolor=OUTLINE, linewidth=0.3))
lax.text(0.05, 0.5, f"{one_gen + one_rna} with one supported call ({one_rna} RNA, {one_gen} genetics)\n"
         f"{both_untestable} untestable in both domains; counted in the RNA-untestable bar",
         ha="left", va="center", linespacing=1.0)
fig.savefig(OUT / "fig6d_next_experiment.pdf")
plt.close(fig)
src6d.append(dict(primary_evidence_class="unresolved", class_display="unresolved", segment="annotation: fine-mapped coding variant PIP >= 0.5",
                  n_genes=n_unres_coding, hatched=False))
src6d.append(dict(primary_evidence_class="all", class_display="total", segment="n genes", n_genes=N, hatched=False))
src6d.append(dict(primary_evidence_class="untested", class_display="untested", segment="annotation: both domains untestable, counted in RNA untestable segment", n_genes=both_untestable, hatched=False))
pd.DataFrame(src6d).to_csv(SRC / "6D.tsv", sep="\t", index=False)
print("6B", m.to_dict(), "coding", mc.to_dict())
print("6D", cc.to_dict(), n_no_rna, n_no_gen, one_gen, one_rna)
