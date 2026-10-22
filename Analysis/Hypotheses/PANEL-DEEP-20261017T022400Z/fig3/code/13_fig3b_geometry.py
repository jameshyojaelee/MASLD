"""Fig 3B: read every printed cell count with its position, map to (NAS, fibrosis) by row/column labels, compare with the
recomputed grid; record fill colour per cell and text colour. Usage: python 13_fig3b_geometry.py <pdf> <tag>"""
import pdfplumber, pandas as pd, numpy as np, sys
pdf,tag=sys.argv[1],sys.argv[2]
p=pdfplumber.open(pdf).pages[0]
w=p.extract_words(extra_attrs=['size','non_stroking_color'])
W=pd.DataFrame(w); W['cx']=(W.x0+W.x1)/2; W['cy']=(W.top+W.bottom)/2
# column labels: fibrosis 0-4 + Total at bottom; row labels NAS 0-8 at left
bottom=W[(W.top>W.top.max()-30)&(W.text.isin(['0','1','2','3','4','Total']))&(W['size']>5.9)]
bottom=bottom[bottom.top==bottom.top[bottom.text=='Total'].max()] if (bottom.text=='Total').any() else bottom
colx={r.text:r.cx for r in bottom.itertuples()}
left=W[(W.x0<20)&(W.text.str.fullmatch(r'\d'))]; rowy={r.text:r.cy for r in left.itertuples()}
top=W[(W.text=='Total')&(W.x0<30)]; 
if len(top): rowy['Total']=top.cy.iloc[0]
print("col x",colx); print("row y",rowy)
g=pd.read_csv("out/fig3b_recomputed_grid.tsv",sep="\t",index_col=0)
g.index=g.index.astype(str); g.columns=[str(c) for c in g.columns]; g=g.rename(index={'All':'Total'},columns={'All':'Total'})
cells=W[W.text.str.fullmatch(r'\d+')&~W.index.isin(left.index)&~W.index.isin(bottom.index)]
rects=[r for r in p.rects if (r['x1']-r['x0'])<20]
out=[]
for r in cells.itertuples():
    col=min(colx,key=lambda k:abs(colx[k]-r.cx)); row=min(rowy,key=lambda k:abs(rowy[k]-r.cy))
    if abs(colx[col]-r.cx)>4 or abs(rowy[row]-r.cy)>4: continue
    exp=g.loc[row,col] if (row in g.index and col in g.columns) else np.nan
    fill=[tuple(round(v,3) for v in (q['non_stroking_color'] if isinstance(q['non_stroking_color'],(list,tuple)) else [q['non_stroking_color']])) for q in rects if q['x0']<=r.cx<=q['x1'] and q['top']<=r.cy<=q['bottom']]
    out.append(dict(nas=row,fib=col,printed=int(r.text),expected=exp,text_color=str(r.non_stroking_color),fill=str(fill[-1]) if fill else ''))
O=pd.DataFrame(out); O['ok']=O.printed==O.expected
print(O.to_string()); print("checked",len(O),"off",(~O.ok).sum())
# empty cells expected zero
for row in g.index:
    for col in g.columns:
        if g.loc[row,col]>0 and not ((O.nas==row)&(O.fib==col)).any(): print("MISSING printed value for",row,col,g.loc[row,col])
O.to_csv(f"out/fig3b_geometry_{tag}.tsv",sep="\t",index=False)
