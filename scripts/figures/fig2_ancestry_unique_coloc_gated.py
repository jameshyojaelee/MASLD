#!/usr/bin/env python
"""fig2_ancestry_unique_coloc_gated.py  (2026-06-25)  -- Fig 2G (GWS-gated variants)

Honest variants of Fig2G that gate the "non-EUR-unique" colocalizing genes on the
genome-wide significance of their colocalizing lead variant in the non-European GWAS
(audit: scripts/figures/audit_noneur_gws.R -> noneur_gws_audit.csv).

The ungated Fig2G (fig2_ancestry_unique_coloc.py) is LEFT UNTOUCHED. This script adds:
  GATE=gws         -> p<5e-8   -> Fig2G_ancestry_unique_coloc_GWS.pdf
  GATE=suggestive  -> p<1e-6   -> Fig2G_ancestry_unique_coloc_suggestive.pdf
Genes that fail the gate stay colocalizing genes but move to a grey
"non-EUR sub-threshold" tile (they are NOT genome-wide-significant discoveries).

ANCESTRY (2026-07-05): now registry-driven (mirrors gwas_ancestry() in
load_figure_data.R), covering the 50-GWAS MVP portfolio incl. AMR.

STALE-AUDIT GUARD (2026-07-05): noneur_gws_audit.csv is produced by the NON-owned,
pre-MVP audit_noneur_gws.R (same retired ancestry heuristic; reads only BBJ/PanUKBB
sumstat paths, so it cannot resolve MVP GWS p-values). Until that audit is regenerated
for the 50-GWAS portfolio it covers only a fraction of the current non-EUR-unique
genes; this script REFUSES to write a gated panel from a stale audit (see guard below).

Run: GATE=gws        ~/micromamba/envs/rnaseq/bin/python scripts/figures/fig2_ancestry_unique_coloc_gated.py
     GATE=suggestive ~/micromamba/envs/rnaseq/bin/python scripts/figures/fig2_ancestry_unique_coloc_gated.py
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
TITLE_FS, BODY_FS, TXT = 10.5, 8.5, "black"

GATE = os.environ.get("GATE", "gws")
THRESH = {"gws": 5e-8, "suggestive": 1e-6}[GATE]
GLAB = {"gws": "GWS-gated (p < 5×10⁻⁸)", "suggestive": "suggestive-gated (p < 10⁻⁶)"}[GATE]
ULAB = {"gws": "non-EUR\nGWS-unique", "suggestive": "non-EUR\nsuggestive-unique"}[GATE]
OUT  = {"gws": "Fig2G_ancestry_unique_coloc_GWS.pdf",
        "suggestive": "Fig2G_ancestry_unique_coloc_suggestive.pdf"}[GATE]

BASE = os.environ.get("MASLD_PROJECT_ROOT",
                      "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
PANEL_DIR = os.path.join(BASE, "figures/main/fig2_genetics/panels")
SC = os.path.join(BASE, "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv")
AUD = os.path.join(BASE, "RNA-seq/results/coloc_variant_classes/noneur_gws_audit.csv")
REG = os.path.join(BASE, "GWAS/finemapping/config/gwas_registry.tsv")

# ancestry from the GWAS registry (mirrors gwas_ancestry() in load_figure_data.R);
# retires the hardcoded name-sets, which had no AMR bin and routed MVP strata into EUR.
_reg = pd.read_csv(REG, sep="\t")
ANC_MAP = dict(zip(_reg["study_name"], _reg["ancestry"]))


def ancestry(g):
    return ANC_MAP.get(g)


# ---- bins (same as the ungated panel) --------------------------------------
df = pd.read_csv(SC)
df["ancestry"] = df["gwas_name"].map(ancestry)
df = df[df["ancestry"].notna()].copy()
df["pp4_best"] = df[["PP.H4.susie", "PP.H4.abf"]].max(axis=1)
sig = df[df["pp4_best"] > 0.5]
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

# STALE-AUDIT GUARD: the audit must cover the CURRENT non-EUR-unique gene set. If it
# was generated on an older (pre-MVP) portfolio it will be missing most genes, and
# every missing gene would be silently dumped into "sub-threshold" -> a wrong panel.
import sys
covered = sum(1 for g in nonEUR_genes if g in aud.index)
frac = covered / max(len(nonEUR_genes), 1)
if frac < 0.9:
    print(f"WARNING/SKIP [{GATE}]: noneur_gws_audit.csv covers only {covered}/"
          f"{len(nonEUR_genes)} ({100*frac:.0f}%) of the current non-EUR-unique genes -- "
          f"it is STALE (pre-MVP). Its generator scripts/figures/audit_noneur_gws.R has "
          f"NOT been re-run for the 50-GWAS MVP portfolio (same retired ancestry heuristic; "
          f"reads only BBJ/PanUKBB sumstat paths, so it cannot resolve MVP GWS p-values). "
          f"Refusing to write {OUT} from a stale audit -- regenerate the audit first.",
          file=sys.stderr)
    sys.exit(0)

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
# light fills so all-black labels stay legible (orange = EAS, matching the ungated panel)
cols  = ["#D6DBDE", "#A7B6BE", "#E9ECEE", "#F79268"]   # grey, grey-blue, light-grey, EAS-orange
W, H = 100.0, 64.0
rects = squarify.squarify(squarify.normalize_sizes(sizes, W, H), 0, 0, W, H)

fig, ax = plt.subplots(figsize=(4.3, 3.0))
labels = [("EUR-only", f"{n_eur} ({pct(n_eur):.0f}%)"),
          ("shared", f"{n_shared} ({pct(n_shared):.0f}%)"),
          ("non-EUR\nsub-threshold", f"{n_fail} ({pct(n_fail):.0f}%)"),
          (ULAB, f"{n_pass} ({pct(n_pass):.0f}%)")]
for i, (r, c, (name, sub)) in enumerate(zip(rects, cols, labels)):
    ax.add_patch(mpatches.Rectangle((r["x"], r["y"]), r["dx"], r["dy"],
                 facecolor=c, edgecolor="white", linewidth=1.4))
    # inside-label the three context tiles only; the discovery tile is small -> leader
    if i < 3 and r["dx"] * r["dy"] >= 32:
        ax.text(r["x"] + r["dx"] / 2, r["y"] + r["dy"] / 2, f"{name}\n{sub}",
                ha="center", va="center", color=TXT, fontsize=BODY_FS, linespacing=1.25)

# always leader the (small) discovery tile to the right margin
rp = rects[3]
ax.annotate(f"{ULAB}\n{n_pass} ({pct(n_pass):.0f}%)", xy=(rp["x"] + rp["dx"] / 2, rp["y"] + rp["dy"] / 2),
            xytext=(W + 3, rp["y"] + rp["dy"] / 2), ha="left", va="center",
            fontsize=BODY_FS, color=TXT,
            arrowprops=dict(arrowstyle="-", lw=0.7, color="black"))
ax.text(W / 2, H + 5.5, f"surviving non-EUR discoveries: {comp_str}",
        ha="center", va="top", fontsize=BODY_FS, color=TXT)

ax.set_xlim(0, W + 30)
ax.set_ylim(0, H + 12)
ax.invert_yaxis()
ax.axis("off")
ax.set_title(f"Colocalization by ancestry — {GLAB}", fontsize=TITLE_FS,
             fontweight="bold", color=TXT, loc="left", pad=6)

os.makedirs(PANEL_DIR, exist_ok=True)
out = os.path.join(PANEL_DIR, OUT)
fig.tight_layout()
fig.savefig(out, bbox_inches="tight")
plt.close(fig)
print(f"[fig2G/{GATE}] wrote {OUT} | EUR-only {n_eur}, shared {n_shared}, "
      f"non-EUR sub-threshold {n_fail}, non-EUR-unique(gated) {n_pass} [{comp_str}]")
print(f"CAPTION: Non-EUR-unique colocalizing genes gated on {GLAB} of the colocalizing lead "
      f"variant in the non-European GWAS (50-GWAS portfolio incl. MVP). Of {len(nonEUR_genes)} "
      f"non-EUR-unique genes, {n_pass} survive ({comp_str}); the remaining {n_fail} are "
      f"sub-threshold (coloc.abf on the eQTL cis-window, not GWS).")
