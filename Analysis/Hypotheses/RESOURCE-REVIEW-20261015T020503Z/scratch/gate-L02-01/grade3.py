import pandas as pd, re, numpy as np
L='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L02-01'
m=pd.read_csv(f'{L}/liver_rna_meta.tsv',sep='\t')
def pct(t):
    vals=[]
    for mm in re.finditer(r'(\d+)\s*(?:-|to)\s*(\d+)\s*\+?\s*%|([<>~]?)\s*(\d+)\s*\+?\s*%',t):
        if mm.group(1): vals.append((float(mm.group(1))+float(mm.group(2)))/2)
        else:
            v=float(mm.group(4)); vals.append(v/2 if mm.group(3)=='<' else v)
    return vals
SEV=r'(mild[- ]to[- ]moderate|mild-moderate|moderate-severe|minimal|mild|moderate|severe|marked|extensive)'
MAP={'minimal':'minimal','mild':'mild','mild to moderate':'mild_moderate','mild-moderate':'mild_moderate','moderate':'moderate','moderate-severe':'severe','severe':'severe','marked':'severe','extensive':'severe'}
ORD={'none':0,'minimal':1,'mild':2,'mild_moderate':2.5,'moderate':3,'severe':4}
rows=[]
for _,r in m.iterrows():
    t=str(r.SMPTHNTS).lower().replace('!0%','10%').replace('%60','60%').replace('80+%','80%').replace('90+ %','90%')
    fat=bool(re.search(r'steato|fat',t))
    neg=bool(re.search(r'no significant steatosis|no steatosis|without steatosis|no fat|no significant fat',t))
    g='no_mention'; p=None
    if fat and neg: g='none'
    elif fat:
        mf=re.search(r'steato|fat',t); i=mf.start()
        win=t[max(0,i-60):min(len(t),i+60)]
        # drop text about fibrous tissue
        pv=[v for v in pct(re.sub(r'fibrous tissue[^;]*','',win))]
        p=max(pv) if pv else None
        if p is not None:
            g='minimal' if p<5 else ('mild' if p<=33 else ('moderate' if p<=66 else 'severe'))
        else:
            before=list(re.finditer(SEV,t[max(0,i-45):i]))
            after=re.search(SEV,t[i:i+30])
            w=before[-1].group(1) if before else (after.group(1) if after else None)
            g=MAP.get(w.replace('mild to moderate','mild to moderate'),'unspecified') if w else 'unspecified'
            if w and w.startswith('mild') and 'moderate' in w: g='mild_moderate'
    fib=bool(re.search(r'fibros|cirrho|steatohepatitis|mallory|hyalin',t))
    rows.append(dict(donor=r.donor,grade=g,ord=ORD.get(g,np.nan),pct=p,fib_or_nash=fib,note=r.SMPTHNTS))
d=pd.DataFrame(rows)
d.to_csv(f'{L}/note_grades.tsv',sep='\t',index=False)
print(d.grade.value_counts().to_string())
print('graded:',d.ord.notna().sum(),' stated pct:',d.pct.notna().sum(),' pct>33:',(d.pct>33).sum())
print('fib/NASH/cirrhosis/Mallory mention: all',d.fib_or_nash.sum(),'; graded>=moderate with it',(d.fib_or_nash&(d.ord>=3)).sum(), '; unspecified with it',(d.fib_or_nash&(d.grade=="unspecified")).sum())
print('UNSPEC',d[d.grade=='unspecified'].note.str.slice(0,100).tolist())
for g in ['moderate','mild']:
    print(g, d[d.grade==g].note.str.slice(0,80).tolist()[:40])
