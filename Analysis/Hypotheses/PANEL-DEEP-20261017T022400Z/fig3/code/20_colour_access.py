"""Colour-vision simulation (Machado, Oliveira & Fernandes 2009, severity 1.0 matrices applied in linear RGB) and WCAG
2.x contrast for the mark classes each Fig 3 panel asks the reader to tell apart. Distances are CIE76 Delta E in CIELAB
(D65). Pairs below 10 after simulation are flagged as hard to tell apart for small marks; non-text graphics against
white below 3:1 are flagged (WCAG 1.4.11); text below 4.5:1 is flagged (WCAG 1.4.3)."""
import numpy as np, itertools, pandas as pd
M={'protan':np.array([[0.152286,1.052583,-0.204868],[0.114503,0.786281,0.099216],[-0.003882,-0.048116,1.051998]]),
   'deutan':np.array([[0.367322,0.860646,-0.227968],[0.280085,0.672501,0.047413],[-0.011820,0.042940,0.968881]]),
   'tritan':np.array([[1.255528,-0.076749,-0.178779],[-0.078411,0.930809,0.147602],[0.004733,0.691367,0.303900]])}
def lin(c): c=np.asarray(c,float); return np.where(c<=0.04045,c/12.92,((c+0.055)/1.055)**2.4)
def delin(c): c=np.clip(c,0,1); return np.where(c<=0.0031308,12.92*c,1.055*c**(1/2.4)-0.055)
def sim(c,k): return delin(M[k]@lin(c)) if k!='normal' else np.asarray(c,float)
def lab(c):
    l=lin(c); X=np.array([[0.4124,0.3576,0.1805],[0.2126,0.7152,0.0722],[0.0193,0.1192,0.9505]])@l/np.array([0.95047,1,1.08883])
    f=np.where(X>0.008856,np.cbrt(X),7.787*X+16/116); return np.array([116*f[1]-16,500*(f[0]-f[1]),200*(f[1]-f[2])])
def lum(c): l=lin(c); return 0.2126*l[0]+0.7152*l[1]+0.0722*l[2]
def cr(a,b): la,lb=sorted([lum(a),lum(b)],reverse=True); return (la+0.05)/(lb+0.05)
h=lambda s:np.array([int(s[i:i+2],16)/255 for i in (1,3,5)])
sets={
 '3A tile states':{'Available':h('#FFCD69'),'Partial':h('#ED87BD'),'Binned':h('#AD93E6'),'Absent':h('#ECECEC')},
 '3A bars':{'Disease':h('#C9265E'),'Control':h('#9E9E9E')},
 '3C study vs class':{'GSE135251':(0.909804,0.741176,0.913725),'GSE130970':(0.788235,0.721569,0.909804),'GSE126848':(0.658824,0.780392,0.898039),
   'GSE162694':(0.952941,0.745098,0.756863),'GSE213621':(0.984314,0.870588,0.741176),'unique':(0.835294,0.835294,0.835294),
   'same2':(0.890196,0.745098,0.909804),'same3':(0.745098,0.67451,0.878431),'same4':(0.639216,0.67451,0.827451),'same5':(0.498039,0.529412,0.713725),'mixed':(0.709804,0.709804,0.709804)},
 '3D fibrosis ramp':{'Control':(0.619608,)*3,'F0':(0.980392,0.815686,0.807843),'F1':(0.929412,0.647059,0.698039),'F2':(0.882353,0.482353,0.588235),'F3':(0.835294,0.313725,0.478431),'F4':(0.788235,0.14902,0.368627)},
 '3E programs':{'P1':(1.0,0.803922,0.411765),'P2':(0.545098,0.756863,0.388235),'P3':(0.164706,0.54902,0.490196),'P4':(0.929412,0.333333,0.396078),'P5':(0.929412,0.529412,0.741176)},
 '3E bars':{'Up':h('#C9265E'),'Down':h('#1565C0')},
 '3G classes':{'Neither':h('#9E9E9E'),'Transcriptomics':h('#C9265E'),'Genetics':h('#1565C0'),'Both':h('#6A3D9A')},
 'lncRNA':{'Up':h('#C9265E'),'Down':h('#1565C0')},
}
rows=[]
for nm,S in sets.items():
    for k in ['normal','deutan','protan','tritan']:
        L={c:lab(sim(v,k)) for c,v in S.items()}
        for a,b in itertools.combinations(S,2): rows.append(dict(set=nm,vision=k,a=a,b=b,dE=np.linalg.norm(L[a]-L[b])))
    for c,v in S.items(): rows.append(dict(set=nm,vision='contrast_vs_white',a=c,b='white',dE=cr(v,(1,1,1))))
D=pd.DataFrame(rows); D.to_csv("out/colour_access.tsv",sep="\t",index=False)
pd.set_option('display.width',200); pd.set_option('display.max_rows',500)
print("pairs with dE<10 (any vision):"); print(D[(D.vision!='contrast_vs_white')&(D.dE<10)].round(1).to_string())
print("non-text contrast vs white < 3:1:"); print(D[(D.vision=='contrast_vs_white')&(D.dE<3)].round(2).to_string())
# 3B text on fills
G=pd.read_csv("out/fig3b_geometry_main.tsv",sep="\t"); print(G.columns.tolist())
import ast
def t(s):
    v=ast.literal_eval(str(s)); v=v if isinstance(v,(tuple,list)) else (v,)
    return np.array(v*3 if len(v)==1 else v,float)
Ob=pd.read_csv("geom/fig3b_nas_fib_grid_objects.tsv",sep="\t"); R=Ob[(Ob.kind=='rect')&(Ob.fill.astype(str)!='0')&((Ob.x1-Ob.x0)<20)]
W=Ob[(Ob.kind=='char')&(Ob.text.astype(str).str.match(r'^[0-9]$'))]
# rebuild cell fill from the last filled rect under each printed number (stroke-only overlays carry no fill)
def cellfill(row):
    m=W[(W.x0>=row.x0)&(W.x1<=row.x1)&(W.top>=row.top)&(W.bottom<=row.bottom)]
    return m
cells=[]
for _,r in R.iterrows():
    m=cellfill(r)
    if len(m): cells.append(dict(text=''.join(m.sort_values('x0').text.astype(str)),text_color=m.fill.iloc[0],fill=r.fill,x0=r.x0,top=r.top))
G=pd.DataFrame(cells).drop_duplicates(['x0','top'],keep='last')
G['cr']=[cr(t(a),t(b)) for a,b in zip(G.text_color,G.fill)]
G['cr_black']=[cr((0,0,0),t(b)) for b in G.fill]
print(G.sort_values('cr').head(10).round(2).to_string())
print("3B cells",len(G),"text below 4.5:1:",(G.cr<4.5).sum(),"min %.2f"%G.cr.min(),"| if all text black: below 4.5:",(G.cr_black<4.5).sum(),"min %.2f"%G.cr_black.min())
G.to_csv("out/fig3b_text_contrast.tsv",sep="\t",index=False)
