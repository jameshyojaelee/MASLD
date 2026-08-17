#!/usr/bin/env python
"""fig2_ancestry_unique_coloc_gated.py  (2026-06-25)  -- Fig 2E (GWS-gated variants)

Honest variants of Fig2G that gate the "non-EUR-unique" colocalizing genes on the
genome-wide significance of their colocalizing lead variant in the non-European GWAS
(audit: scripts/figures/audit_noneur_gws.R -> noneur_gws_audit.csv).

COLOC SET: multi-signal COLOC. Colocalization is defined as PP.H4.susie > 0.5
ONLY (462 named genes in the promoted 2026-08-17 release), NOT the SuSiE-OR-ABF union that this script used
before. The ancestry partition (EUR-only / shared / non-EUR-unique) is recomputed on the
SuSiE set, so the counts here are SMALLER than the retired union panel.

The ungated Fig2G base panel (fig2_ancestry_unique_coloc.py) is RETIRED 2026-07-06
(indefensible without a non-EUR significance gate); this script produces the gated
SuSiE-COLOC companions that supersede it:
  GATE=gws         -> p<5e-8   -> Fig2E_ancestry_unique_coloc_GWS.pdf
  GATE=suggestive  -> p<1e-6   -> Fig2E_ancestry_unique_coloc_suggestive.pdf
Genes that fail the gate stay SuSiE-coloc genes but move to a grey
"non-EUR sub-threshold" tile (they are NOT genome-wide-significant discoveries).

ANCESTRY (2026-07-05): now registry-driven (mirrors gwas_ancestry() in
load_figure_data.R), covering the 50-GWAS MVP portfolio incl. AMR.

RELEASE GUARD (2026-08-17): the GWS audit gene set must exactly match the
non-EUR-unique set derived from the supplied promoted COLOC table. Coverage-only
agreement is insufficient because membership changed at promotion.

Run: GATE=gws        ~/micromamba/envs/rnaseq/bin/python scripts/figures/fig2_ancestry_unique_coloc_gated.py
     GATE=suggestive ~/micromamba/envs/rnaseq/bin/python scripts/figures/fig2_ancestry_unique_coloc_gated.py
"""
import os
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import squarify
import pandas as pd

plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Nimbus Sans", "DejaVu Sans"],
                     "pdf.fonttype": 42, "axes.linewidth": 0})
# consistent across all 3 ancestry-unique figures: Helvetica, all-black text.
# ALL text = 6 pt, non-bold (lab style, 2026-07-06): title AND every tile/caption/number
# label share one size so nothing outsizes the title (matches the other Fig 2 panels).
TITLE_FS, BODY_FS, TXT = 6, 6, "black"

GATE = os.environ.get("GATE", "gws")
THRESH = {"gws": 5e-8, "suggestive": 1e-6}[GATE]
GLAB = {"gws": "GWS-gated (p < 5×10⁻⁸)", "suggestive": "suggestive-gated (p < 10⁻⁶)"}[GATE]
ULAB = {"gws": "non-EUR\nGWS-unique", "suggestive": "non-EUR\nsuggestive-unique"}[GATE]
OUT  = {"gws": "Fig2E_ancestry_unique_coloc_GWS.pdf",
        "suggestive": "Fig2E_ancestry_unique_coloc_suggestive.pdf"}[GATE]

BASE = os.environ.get("MASLD_PROJECT_ROOT",
                      "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
PANEL_DIR = os.environ.get(
    "FIG2_CANDIDATE_DIR", os.path.join(BASE, "figures/main/fig2_genetics/panels")
)
SC = os.environ.get(
    "FIG2_COLOC_INPUT",
    os.path.join(BASE, "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv"),
)
AUD = os.environ.get(
    "FIG2_NONEUR_AUDIT",
    os.path.join(BASE, "RNA-seq/results/coloc_variant_classes/noneur_gws_audit.csv"),
)
REG = os.path.join(BASE, "GWAS/finemapping/config/gwas_registry.tsv")
TIER = os.path.join(BASE, "GWAS/finemapping/config/gwas_trait_tier.tsv")

# ancestry from the GWAS registry (mirrors gwas_ancestry() in load_figure_data.R);
# retires the hardcoded name-sets, which had no AMR bin and routed MVP strata into EUR.
_reg = pd.read_csv(REG, sep="\t")
ANC_MAP = dict(zip(_reg["study_name"], _reg["ancestry"]))

# MAIN (Tier-1/2, liver-specific) allowlist (2026-07-06): keep only placement=="main"
# strata so the gated non-EUR-unique set MATCHES the MAIN-scoped ungated Fig2G panel and
# the MAIN-scoped audit (audit_noneur_gws.R). Tier-3/4 supp strata are excluded.
_tier = pd.read_csv(TIER, sep="\t")
MAIN_STUDIES = set(_tier.loc[_tier["placement"] == "main", "study_name"])


def ancestry(g):
    return ANC_MAP.get(g)


# ---- bins (SuSiE-COLOC set; partition mirrors the retired ungated panel) ----
df = pd.read_csv(SC)
df = df[df["gwas_name"].isin(MAIN_STUDIES)].copy()   # MAIN (Tier-1/2) strata only
df["ancestry"] = df["gwas_name"].map(ancestry)
df = df[df["ancestry"].notna()].copy()
# Multi-signal COLOC: PP.H4.susie > 0.5 in the promoted input, not the
# historical 473-gene July set. NaN PP.H4.susie is excluded by the comparison.
df["pp4_best"] = df["PP.H4.susie"]
sig = df[df["pp4_best"] > 0.5]
sig = sig[sig["gene"].notna() & sig["gene"].astype(str).str.strip().ne("")]
by_gene = sig.groupby("gene")["ancestry"].apply(set)

def binof(a):
    nonE = a - {"EUR"}
    if a == {"EUR"}: return "EUR only"
    if "EUR" in a and nonE: return "shared"
    if "EUR" not in a: return "non-EUR-unique"
    return None
gb = by_gene.apply(binof)
n_eur = int((gb == "EUR only").sum())
n_shared = int((gb == "shared").sum())
nonEUR_genes = set(gb[gb == "non-EUR-unique"].index)
n_total = int(by_gene.shape[0])

# ---- apply the GWS / suggestive gate to the non-EUR-unique genes ------------
aud = pd.read_csv(AUD).set_index("gene")

# Release guard: exact gene-set agreement is required, not merely high coverage.
import sys
audit_genes = set(aud.index.astype(str))
missing = sorted(nonEUR_genes - audit_genes)
extra = sorted(audit_genes - nonEUR_genes)
if missing or extra:
    print(f"ERROR [{GATE}]: non-EUR GWS audit does not exactly match the promoted "
          f"non-EUR-unique set (missing={len(missing)}, extra={len(extra)}). "
          "Regenerate it with audit_noneur_gws.R against the same FIG2_COLOC_INPUT. "
          f"Missing examples: {missing[:5]}; extra examples: {extra[:5]}",
          file=sys.stderr)
    sys.exit(1)

minp = aud["min_p"].to_dict()
pass_genes = [g for g in nonEUR_genes if pd.notna(minp.get(g)) and minp[g] < THRESH]
fail_genes = [g for g in nonEUR_genes if g not in pass_genes]
n_pass, n_fail = len(pass_genes), len(fail_genes)
# ancestry composition of the surviving discoveries (from the audit bin)
comp = aud.loc[[g for g in pass_genes if g in aud.index], "bin"].value_counts().to_dict()
comp_str = ", ".join(f"{v} {k.replace(' only','')}" for k, v in sorted(comp.items(), key=lambda x: -x[1]))
pct = lambda n: 100.0 * n / n_total

# ---- flat treemap: EUR-only / shared / non-EUR sub-threshold / non-EUR-unique
sizes = [n_eur, n_shared, n_fail, n_pass]
# light fills so all-black labels stay legible (orange = surviving non-EUR discovery tile)
cols  = ["#D6DBDE", "#A7B6BE", "#E9ECEE", "#F79268"]   # grey, grey-blue, light-grey, EAS-orange
W, H = 100.0, 64.0
rects = squarify.squarify(squarify.normalize_sizes(sizes, W, H), 0, 0, W, H)

fig, ax = plt.subplots(figsize=(2.00, 2.10))
labels = [("EUR-only", f"{n_eur} ({pct(n_eur):.0f}%)"),
          ("shared", f"{n_shared} ({pct(n_shared):.0f}%)"),
          ("non-EUR\nsub-threshold", f"{n_fail} ({pct(n_fail):.0f}%)"),
          (ULAB, f"{n_pass} ({pct(n_pass):.0f}%)")]
for i, (r, c, (name, sub)) in enumerate(zip(rects, cols, labels)):
    ax.add_patch(mpatches.Rectangle((r["x"], r["y"]), r["dx"], r["dy"],
                 facecolor=c, edgecolor="white", linewidth=1.4))
    # inside-label the three context tiles only; the discovery tile is small -> leader
    if i < 2 and r["dx"] * r["dy"] >= 32:
        ax.text(r["x"] + r["dx"] / 2, r["y"] + r["dy"] / 2, f"{name}\n{sub}",
                ha="center", va="center", color=TXT, fontsize=BODY_FS, linespacing=1.25)

# Leader the two small tiles into separate whitespace regions below the mosaic.
rf = rects[2]
ax.annotate(f"non-EUR sub-threshold\n{n_fail} ({pct(n_fail):.0f}%)",
            xy=(rf["x"] + rf["dx"] / 2, rf["y"] + rf["dy"] / 2),
            xytext=(27, H + 4.0), ha="center", va="top",
            fontsize=BODY_FS, color=TXT,
            arrowprops=dict(arrowstyle="-", lw=0.45, color="black",
                            shrinkA=0, shrinkB=0))
rp = rects[3]
ax.annotate(f"{ULAB}\n{n_pass} ({pct(n_pass):.0f}%)", xy=(rp["x"] + rp["dx"] / 2, rp["y"] + rp["dy"] / 2),
            xytext=(75, H + 4.0), ha="center", va="top",
            fontsize=BODY_FS, color=TXT,
            arrowprops=dict(arrowstyle="-", lw=0.45, color="black",
                            shrinkA=0, shrinkB=0))
ax.text(W / 2, H + 13.0, f"multi-signal COLOC\nnon-EUR GWS: {comp_str}",
        ha="center", va="top", fontsize=BODY_FS, linespacing=1.1, color=TXT)

ax.set_xlim(0, W)
ax.set_ylim(0, H + 19)
ax.invert_yaxis()
ax.axis("off")

os.makedirs(PANEL_DIR, exist_ok=True)
out = os.path.join(PANEL_DIR, OUT)
fig.subplots_adjust(left=0.04, right=0.96, top=0.97, bottom=0.03)
fig.savefig(out)   # NO bbox_inches="tight": exact figsize for place-at-100%
plt.close(fig)
source_out = os.path.join(PANEL_DIR, "Fig2E_ancestry_unique_coloc_GWS_source.tsv")
with open(source_out, "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["category", "n_genes", "percent_of_multi_signal_coloc", "threshold"])
    for category, count in zip(
        ["EUR_only", "shared", "non_EUR_sub_threshold", "non_EUR_GWS_unique"], sizes
    ):
        writer.writerow([category, count, f"{pct(count):.6f}",
                         "SuSiE PP.H4 > 0.5; non-EUR GWS p < 5e-8"])
print(f"[fig2E/{GATE}] wrote {OUT} | multi-signal COLOC (SuSiE PP.H4>0.5): EUR-only {n_eur}, "
      f"shared {n_shared}, non-EUR sub-threshold {n_fail}, non-EUR-unique(gated) {n_pass} [{comp_str}]")
print(f"CAPTION: Non-EUR-unique multi-signal COLOC genes (SuSiE PP.H4 > 0.5) gated on {GLAB} of the "
      f"colocalizing lead variant in the non-European GWAS (35 Tier-1/2 liver-specific GWAS incl. "
      f"MVP NAFLD/ALT/AST). Of {len(nonEUR_genes)} non-EUR-unique genes, {n_pass} survive "
      f"({comp_str}); the remaining {n_fail} are sub-threshold (multi-signal COLOC but the non-EUR lead "
      f"variant is not genome-wide significant).")
