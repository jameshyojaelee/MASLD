"""Fig 4D geometry/label check: read each cell's printed count and fill; compare to r4d_counts.tsv.
Usage: python g4d_geometry.py <pdf> <tag>"""
import sys
from pdfmarks import *
pdfp, tag = sys.argv[1], sys.argv[2]
out = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig4/recompute"
X = pd.read_csv(f"{out}/r4d_counts.tsv", sep="\t", index_col=0)
pdf, pg, ch, w, r, l, c = load(pdfp)
rows = ["Hep", "Fib", "Chol", "Mac", "T"]; cols = ["Hep", "Fib", "Chol", "Mac", "T", "Endo"]
# row labels: words left of the grid; column labels: words below grid
num = w[w.text.str.fullmatch(r"\d+")]
labs = w[w.text.isin(rows + cols)]
rowlab = labs[labs.x1 < num.x0.min()].copy()
collab = labs[labs.top > num.bottom.max()].copy()
res = []
rowy = {t.text: t.yc for _, t in rowlab.iterrows()}
colx = {t.text: t.xc for _, t in collab.iterrows()}
print("row labels", rowy, "\ncol labels", colx)
cells = r[(r.w > 10) & (r.w < 60) & (r.h > 5) & (r.h < 30)]
for rn in rows:
    for cn in cols:
        if rn not in rowy or cn not in colx:
            res.append(dict(tag=tag, mark="cell", row=rn, col=cn, expected=int(X.loc[rn, cn]), drawn="label missing", off=True)); continue
        n = num[((num.yc - rowy[rn]).abs() < 4) & ((num.xc - colx[cn]).abs() < 8)]
        val = int(n.text.iloc[0]) if len(n) == 1 else None
        tcol = ch[(ch.text.str.isdigit()) & ((((ch.top + ch.bottom) / 2) - rowy[rn]).abs() < 4) & ((((ch.x0 + ch.x1) / 2) - colx[cn]).abs() < 8)].fill.unique().tolist()
        cc = cells[(cells.x0 < colx[cn]) & (cells.x1 > colx[cn]) & (cells.top < rowy[rn]) & (cells.bottom > rowy[rn])]
        fill = cc.fill.iloc[-1] if len(cc) else "none"
        exp = int(X.loc[rn, cn])
        res.append(dict(tag=tag, mark="cell", row=rn, col=cn, expected=exp, drawn=val, off=(val != exp), fill=fill, text_colour=",".join(tcol),
                        diag=(rn == cn)))
R = pd.DataFrame(res)
R.to_csv(f"{out}/g4d_{tag}.tsv", sep="\t", index=False)
print(R.pivot(index="row", columns="col", values="drawn").reindex(index=rows, columns=cols))
print("cells checked", len(R), "off", int(R.off.sum()))
print(R[R.off].to_string())
print(R.groupby(["diag", "fill", "text_colour"]).size())
print("diag sum", R[R.diag].drawn.sum(), "off-diag sum", R[~R.diag].drawn.sum())
# axis titles / n presence
txt = " ".join(w.text)
rot = "".join(ch[~ch.upright.astype(bool)].sort_values("top", ascending=False).text)
print("rotated text (bottom to top):", rot, "| sorted letters match 'Discovery lineage':", sorted(rot.replace(" ", "")) == sorted("Discoverylineage"))
print("axis titles present: row title", ("Discovery" in txt) or sorted(rot.replace(" ", "")) == sorted("Discoverylineage"), "| column title", "Content-associated" in txt, "; n=117 printed:", "117" in txt)
