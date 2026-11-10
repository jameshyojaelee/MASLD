"""Fig 2E candidate redraw (decision-dependent: draws 'no non-EUR SuSiE result' as its own state).

Edited from code/originals/fig2_ancestry_unique_coloc_gated.py. Kept: the 436 multi-signal genes
(SuSiE PP.H4 > 0.5, 35 tier-1/2 studies) split by the ancestries in which they are supported, the
grey / grey-blue / orange fills, 6 pt black text.
Changed (deep review 2026-10-17, PANEL_SPECS 2E):
  - treemap replaced by horizontal bars on a gene-count axis
  - EUR only split into 'SuSiE run in non-EUR, PP.H4 <= 0.5' (127) and 'no non-EUR SuSiE result' (213),
    the second drawn as its own state, not as negative; genes with non-EUR ABF PP.H4 > 0.5 marked
  - non-EUR only split into lead p < 5e-8 without a EUR multi-signal gene within 1 Mb (5), with one (4),
    and lead p >= 5e-8 (23); GWS defined on the panel; the word 'unique' dropped
Reads source_tables/2E.tsv only.
"""
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd

OUT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CATALOG-ADOPTION-20261107T022216Z/figures"
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Nimbus Sans", "DejaVu Sans"],
    "font.size": 6, "pdf.fonttype": 42, "axes.linewidth": 0.5, "xtick.major.width": 0.5,
    "xtick.major.size": 2, "ytick.major.size": 0, "hatch.linewidth": 0.4, "text.color": "black",
    "xtick.color": "black", "ytick.color": "black", "axes.labelcolor": "black",
})
e = pd.read_csv(f"{OUT}/source_tables/2E.tsv", sep="\t")
STYLE = {  # segment -> (fill, hatch, edge, key label)
    "tested in non-EUR, not supported": ("#A7B6BE", "", "black", "finite non-EUR SuSiE result, PP.H4 ≤ 0.5"),
    "not tested in non-EUR": ("#FFFFFF", "////////", "black", "no finite non-EUR SuSiE result"),
    "shared": ("#D6DBDE", "", "black", "supported in EUR and non-EUR"),
    "lead p < 5e-8, no EUR gene within 1 Mb": ("#F79268", "", "black",
                                               "lead p < 5e-8, no EUR multi-signal gene within 1 Mb"),
    "lead p < 5e-8, EUR gene within 1 Mb": ("#F79268", "......", "black",
                                            "lead p < 5e-8, EUR multi-signal gene within 1 Mb"),
    "lead p >= 5e-8": ("#E9ECEE", "", "black", "lead p ≥ 5e-8"),
}
GROUPS = ["EUR only", "Shared EUR and non-EUR", "Non-EUR only"]
YLAB = {"EUR only": "EUR only", "Shared EUR and non-EUR": "EUR and\nnon-EUR", "Non-EUR only": "non-EUR\nonly"}
H = 0.62
fig = plt.figure(figsize=(2.60, 2.55))
ax = fig.add_axes([0.20, 0.52, 0.74, 0.36])
for i, g in enumerate(GROUPS):
    left = 0
    for _, r in e[e.group == g].iterrows():
        fc, hat, ec, _ = STYLE[r.segment]
        ax.barh(i, r.n_genes, left=left, height=H, color=fc, hatch=hat, edgecolor=ec, linewidth=0.3)
        k = r["of_which_nonEUR_ABF_gt_0.5"]
        if pd.notna(k) and k > 0:  # genes with non-EUR ABF PP.H4 > 0.5: black strip along the bar base
            ax.barh(i + H / 2 - 0.07, k, left=left, height=0.14, color="black", linewidth=0)
        left += r.n_genes
    ax.text(left + 6, i, f"{int(left)}", va="center", ha="left")
ax.set_yticks(range(3))
ax.set_yticklabels([YLAB[g] for g in GROUPS])
ax.set_ylim(2.5, -0.5)
ax.set_xlim(0, 380)
ax.set_xticks([0, 100, 200, 300])
ax.set_xlabel("genes")
for s in ["top", "right"]:
    ax.spines[s].set_visible(False)
ax.set_title("Multi-signal COLOC genes (SuSiE PP.H4 > 0.5),\nn = 436, by ancestry of the supporting GWAS",
             loc="left", fontsize=6, pad=3)
handles = []
for seg in STYLE:
    n = int(e.loc[e.segment == seg, "n_genes"].iloc[0])
    fc, hat, ec, lab = STYLE[seg]
    handles.append(mpatches.Patch(facecolor=fc, hatch=hat, edgecolor=ec, linewidth=0.3, label=f"{lab} ({n})"))
nabf = e["of_which_nonEUR_ABF_gt_0.5"].fillna(0).astype(int)
handles.append(mpatches.Patch(facecolor="black", linewidth=0,
                              label=f"of EUR only: non-EUR ABF PP.H4 > 0.5 ({nabf.iloc[0]} + {nabf.iloc[1]})"))
fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.01, 0.0), ncol=1, frameon=False, fontsize=6,
           handlelength=1.2, handleheight=0.8, handletextpad=0.4, labelspacing=0.25, borderaxespad=0.2,
           title="lead p = GWAS p at the COLOC lead variant\nin the non-EUR study", title_fontsize=6,
           alignment="left")
for txt in fig.findobj(matplotlib.text.Text):
    txt.set_color("black")
out = f"{OUT}/panels/Fig2E_ancestry_unique_coloc_GWS.pdf"
fig.savefig(out, facecolor="white")
print("Saved", out)
