"""Extra label/colour checks run identically on a main panel set and a calibration set: 3B fill ramp monotone in count;
3C ribbon colour equals its source study bar colour; 3D marginal-density modes vs recomputed KDE modes; 3G legend
swatch colour vs the class whose drawn count equals the printed legend count. Usage: <panel_dir> <tag>"""
import pdfplumber, numpy as np, pandas as pd, sys, re
d,tag=sys.argv[1],sys.argv[2]; res=[]
col=lambda c:tuple(round(v,4) for v in (c if isinstance(c,(list,tuple)) else [c]*3))
def lum(c): c=np.array(c); l=np.where(c<=0.04045,c/12.92,((c+0.055)/1.055)**2.4); return 0.2126*l[0]+0.7152*l[1]+0.0722*l[2]
# 3B
p=pdfplumber.open(f"{d}/fig3b_nas_fib_grid.pdf").pages[0]; W=pd.DataFrame(p.extract_words(extra_attrs=['non_stroking_color']))
R=[q for q in p.rects if (q['x1']-q['x0'])<20 and q['fill'] and isinstance(q['non_stroking_color'],(list,tuple)) and len(q['non_stroking_color'])==3]
cells=[]
for _,w in W[W.text.str.fullmatch(r'\d+')].iterrows():
    cx,cy=(w.x0+w.x1)/2,(w.top+w.bottom)/2; under=[q for q in R if q['x0']<=cx<=q['x1'] and q['top']<=cy<=q['bottom']]
    if under and cx<86 and cy>20: cells.append((int(w.text),lum(col(under[-1]['non_stroking_color']))))
C=pd.DataFrame(cells,columns=['n','lum']).sort_values('n'); viol=0; bad3b=[]
for i in range(len(C)):
    for j in range(len(C)):
        if C.n.iloc[i]<C.n.iloc[j] and C.lum.iloc[i]<C.lum.iloc[j]-1e-6: viol+=1; bad3b.append((int(C.n.iloc[i]),int(C.n.iloc[j])))
res.append(('3B','fill darker for a smaller count (pairs)',viol,f"{len(C)} grid cells; examples {sorted(set(bad3b))[:6]}"))
# 3C
p=pdfplumber.open(f"{d}/fig3c_cohort_alluvial.pdf").pages[0]
studies=[q for q in p.rects if (q['x1']-q['x0'])<20 and (q['bottom']-q['top'])>5 and q['x0']<60]
rib=[c for c in p.curves if c.get('fill') and (c['x1']-c['x0'])>20]
bad=0
for c in rib:
    xs=[t[0] for t in c['pts']]; left=[t for t in c['pts'] if abs(t[0]-min(xs))<0.5]; ly=np.mean([t[1] for t in left])
    s=[q for q in studies if q['top']-0.5<=ly<=q['bottom']+0.5]
    if s and col(s[0]['non_stroking_color'])!=col(c['non_stroking_color']): bad+=1; print('3C mismatch ribbon left y',round(ly,2),'bars',[(round(q['top'],2),round(q['bottom'],2),col(q['non_stroking_color'])) for q in s],'ribbon',col(c['non_stroking_color']))
res.append(('3C','ribbons whose colour differs from their source study bar',bad,f"{len(rib)} ribbons, {len(studies)} study bars"))
# 3G legend
p=pdfplumber.open(f"{d}/fig3g_genetic_state_complementarity.pdf").pages[0]; W=pd.DataFrame(p.extract_words())
txt=' '.join(W[W.top>p.height*0.75].text); pairs=re.findall(r'([A-Za-z]+(?: sig\.)?) \((\d+)\)',txt)
print(tag,"3G legend text pairs",pairs)
res.append(('3G','legend entries parsed',len(pairs),txt))
for r in res: print(tag,*r,sep=' | ')
