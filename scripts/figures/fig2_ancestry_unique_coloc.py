#!/usr/bin/env python
"""fig2_ancestry_unique_coloc.py  (2026-06-25)  -- Fig 2G

Ancestry SPECIFICITY of colocalizing genes as a NESTED TREEMAP. Shows EVERYTHING
including EUR (EUR-dominance is fine -> it is honestly the largest tile); area
encodes gene count so the panel fills completely (no wasted whitespace) and the
small non-EUR categories still get a labeled tile instead of an invisible sliver.

  Top level : EUR-only / EUR+non-EUR (shared) / non-EUR-unique
  The non-EUR-unique tile is subdivided by ancestry: EAS / AFR / SAS / AMR / multiple non-EUR
  (35 Tier-1/2 liver-specific GWAS incl. MVP NAFLD/ALT/AST; scoped to MAIN strata 2026-07-06).

Purely GWAS x eQTL colocalization (no RNA-seq disease signal) -> safe before Fig 3.
Colocalization = best PP.H4 > 0.5 (SuSiE with ABF fallback, since SuSiE cannot
converge on the sparse non-EUR PanUKBB/BBJ GWAS). Mirrors the squarify/matplotlib
treemap approach of scripts/figures/dataset_treemap.py.

Run: ~/micromamba/envs/rnaseq/bin/python scripts/figures/fig2_ancestry_unique_coloc.py
Out: figures/main/fig2_genetics/panels/Fig2G_ancestry_unique_coloc.pdf (+ source CSV)
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import squarify
import pandas as pd

plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Nimbus Sans", "DejaVu Sans"],
                     "pdf.fonttype": 42, "axes.linewidth": 0})
# consistent across all 3 ancestry-unique figures: Helvetica, all-black text
TITLE_FS, BODY_FS, TXT = 6, 6, "black"

BASE = os.environ.get("MASLD_PROJECT_ROOT",
                      "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
PANEL_DIR = os.path.join(BASE, "figures/main/fig2_genetics/panels")
SC = os.path.join(BASE, "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv")
REG = os.path.join(BASE, "GWAS/finemapping/config/gwas_registry.tsv")
TIER = os.path.join(BASE, "GWAS/finemapping/config/gwas_trait_tier.tsv")

# ancestry from the GWAS registry (SINGLE SOURCE OF TRUTH). The .py panels cannot
# source load_figure_data.R, so mirror its registry-driven gwas_ancestry() here.
# Retires the hardcoded EUR/EAS/AFR/SAS name-sets, which had NO AMR bin and silently
# routed every MVP stratum (MVP_*_AMR/AFR/EAS/EUR) into "EUR" once the portfolio grew
# 23 -> 50 GWAS with MVP. Ancestry hexes below mirror ANCESTRY_COLORS in that file.
_reg = pd.read_csv(REG, sep="\t")
ANC_MAP = dict(zip(_reg["study_name"], _reg["ancestry"]))

# MAIN (Tier-1/2, liver-specific) allowlist (2026-07-06): keep only placement=="main"
# strata (NAFLD/NASH/PDFF + ALT/AST/GGT); the Tier-3/4 supp strata (MVP Cirrhosis/
# ChronLiver/Albumin/Platelet) move to a supplementary full-portfolio figure.
_tier = pd.read_csv(TIER, sep="\t")
MAIN_STUDIES = set(_tier.loc[_tier["placement"] == "main", "study_name"])


def ancestry(g):
    return ANC_MAP.get(g)


# ---- data: per-gene ancestry set -> specificity bin ------------------------
df = pd.read_csv(SC)
df = df[df["gwas_name"].isin(MAIN_STUDIES)].copy()   # MAIN (Tier-1/2) strata only
df["ancestry"] = df["gwas_name"].map(ancestry)
df = df[df["ancestry"].notna()].copy()
df["pp4_best"] = df[["PP.H4.susie", "PP.H4.abf"]].max(axis=1)
sig = df[df["pp4_best"] > 0.5]

by_gene = sig.groupby("gene")["ancestry"].apply(lambda s: set(s))
def binof(ancs):
    nonE = ancs - {"EUR"}
    if ancs == {"EUR"}:            return "EUR only"
    if "EUR" in ancs and nonE:     return "EUR + non-EUR (shared)"
    if ancs == {"EAS"}:            return "EAS only"
    if ancs == {"AFR"}:            return "AFR only"
    if ancs == {"SAS"}:            return "SAS only"
    if ancs == {"AMR"}:            return "AMR only"
    if "EUR" not in ancs and len(nonE) >= 2: return "multiple non-EUR"
    return None
bins = by_gene.apply(binof).value_counts().to_dict()

n_eur    = bins.get("EUR only", 0)
n_shared = bins.get("EUR + non-EUR (shared)", 0)
eas, afr, sas, amr, multi = (bins.get("EAS only", 0), bins.get("AFR only", 0),
                             bins.get("SAS only", 0), bins.get("AMR only", 0),
                             bins.get("multiple non-EUR", 0))
n_unique = eas + afr + sas + amr + multi
n_total  = int(by_gene.shape[0])
pct = lambda n: 100.0 * n / n_total

# ---- nested squarify layout ------------------------------------------------
W, H = 100.0, 64.0
top_sizes = [n_eur, n_shared, n_unique]                       # EUR-only, shared, non-EUR-unique
top = squarify.squarify(squarify.normalize_sizes(top_sizes, W, H), 0, 0, W, H)
r_eur, r_shared, r_non = top[0], top[1], top[2]
sub_sizes = [eas, afr, sas, amr, multi]
sub = squarify.squarify(
    squarify.normalize_sizes(sub_sizes, r_non["dx"], r_non["dy"]),
    r_non["x"], r_non["y"], r_non["dx"], r_non["dy"])

# ancestry hexes mirror ANCESTRY_COLORS (load_figure_data.R): EAS red / AFR green /
# SAS purple / AMR orange; multi = neutral grey. EUR/shared tiles stay neutral grey.
GREY_E, GREY_S = "#D6DBDE", "#A7B6BE"
ANC = {"EAS": "#C44E52", "AFR": "#55A868", "SAS": "#8172B3",
       "AMR": "#DD8452", "multi": "#9AA7B1"}

fig, ax = plt.subplots(figsize=(4.3, 3.0))


def tile(r, color, lines, fs=BODY_FS):
    ax.add_patch(mpatches.Rectangle((r["x"], r["y"]), r["dx"], r["dy"],
                 facecolor=color, edgecolor="white", linewidth=1.4))
    if r["dx"] * r["dy"] >= 28:                       # only label if the tile fits text
        ax.text(r["x"] + r["dx"] / 2, r["y"] + r["dy"] / 2, "\n".join(lines),
                ha="center", va="center", color=TXT, fontsize=fs, linespacing=1.25)


tile(r_eur,    GREY_E, ["EUR-only", f"{n_eur} ({pct(n_eur):.0f}%)"])
tile(r_shared, GREY_S, ["shared", f"{n_shared} ({pct(n_shared):.0f}%)"])
sub_meta = [("EAS", eas, "EAS"), ("AFR", afr, "AFR"), ("SAS", sas, "SAS"),
            ("AMR", amr, "AMR"), ("multi", multi, "multi")]
for r, (name, n, key) in zip(sub, sub_meta):
    tile(r, ANC[key], [name, str(n)])

# outline groups the non-EUR-unique tiles; header sits in the bottom margin
ax.add_patch(mpatches.Rectangle((r_non["x"], r_non["y"]), r_non["dx"], r_non["dy"],
             facecolor="none", edgecolor="black", linewidth=1.8))
ax.text(r_non["x"] + r_non["dx"] / 2, H + 5.5,
        f"non-EUR-unique\n{n_unique} ({pct(n_unique):.0f}%)",
        ha="center", va="top", color=TXT, fontsize=BODY_FS, linespacing=1.2)
# external leader for sub-tiles too small to hold their own label (e.g. multi)
rx = r_non["x"] + r_non["dx"]
for r, (name, n, key) in zip(sub, sub_meta):
    if r["dx"] * r["dy"] < 28:
        ax.annotate(f"{name} {n}", xy=(r["x"] + r["dx"] / 2, r["y"] + r["dy"] / 2),
                    xytext=(rx + 2.5, r["y"] + r["dy"] / 2),
                    ha="left", va="center", fontsize=BODY_FS, color=TXT,
                    arrowprops=dict(arrowstyle="-", lw=0.6, color="black"))

ax.set_xlim(0, W + 22)
ax.set_ylim(0, H + 13)
ax.invert_yaxis()
ax.axis("off")
print("[caption] Colocalization by ancestry (SuSiE or ABF)")

os.makedirs(PANEL_DIR, exist_ok=True)
out = os.path.join(PANEL_DIR, "Fig2G_ancestry_unique_coloc.pdf")
fig.tight_layout()
# ---- base ungated Fig2G RETIRED 2026-07-06 ---------------------------------
# The ungated base panel counts ALL non-EUR-unique coloc genes with NO non-EUR
# significance gate, which is indefensible (a coloc can score PP.H4 > 0.5 on a
# sub-threshold non-EUR signal). Its output is DISABLED so the retired PDF is not
# regenerated; the retired PDF is archived at
#   figures/main/fig2_genetics/panels/_archive/RETIRED_Fig2G_ancestry_unique_coloc_base.pdf
# Use the GWS/suggestive gated SuSiE-COLOC panels instead
#   (scripts/figures/fig2_ancestry_unique_coloc_gated.py).
# fig.savefig(out, bbox_inches="tight")
# pd.DataFrame({"bin": ["EUR only", "EUR + non-EUR (shared)", "EAS only", "AFR only",
#                       "SAS only", "AMR only", "multiple non-EUR"],
#               "n": [n_eur, n_shared, eas, afr, sas, amr, multi]}).to_csv(
#     os.path.join(PANEL_DIR, "Fig2G_ancestry_unique_coloc_source.csv"), index=False)
plt.close(fig)
print("[fig2G] BASE UNGATED PANEL RETIRED 2026-07-06 — output disabled; "
      "use the GWS/suggestive gated SuSiE-COLOC panels "
      "(fig2_ancestry_unique_coloc_gated.py) instead.")
