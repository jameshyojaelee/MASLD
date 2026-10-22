import pandas as pd, re, numpy as np
L='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L02-01'
m=pd.read_csv(f'{L}/liver_rna_meta.tsv',sep='\t')
def pct(t):
    vals=[]
    for mm in re.finditer(r'(\d+)\s*(?:-|to)\s*(\d+)\s*\+?\s*%|([<>~]?)\s*(\d+)\s*\+?\s*%',t):
        if mm.group(1): vals.append((float(mm.group(1))+float(mm.group(2)))/2)
        else:
            v=float(mm.group(4))
            if mm.group(3)=='<': v=v/2
            vals.append(v)
    return vals
ORD={'none':0,'minimal':1,'mild':2,'mild_moderate':2.5,'moderate':3,'severe':4}
rows=[]
for _,r in m.iterrows():
    t=str(r.SMPTHNTS).lower().replace('!0%','10%').replace('%60','60%').replace('80+%','80%').replace('90+ %','90%')
    fat=bool(re.search(r'steato|fat',t))
    neg=bool(re.search(r'no significant steatosis|no steatosis|without steatosis|no fat|no significant fat',t))
    segs=[s for s in re.split(r'[;,.]',t) if re.search(r'steato|fat',s)]
    st=' '.join(segs)
    pv=[]
    for s in segs: pv+=pct(s)
    if not pv and re.search(r'steato|fat',t) and re.search(r'\d\s*%',t):
        pass
    p=max(pv) if pv else None
    if not fat: g='no_mention'
    elif neg: g='none'
    elif p is not None:
        g='minimal' if p<5 else ('mild' if p<=33 else ('moderate' if p<=66 else 'severe'))
    elif re.search(r'severe|marked|extensive|panlobular',st): g='severe'
    elif re.search(r'mild[- ]to[- ]moderate|mild-moderate',st): g='mild_moderate'
    elif re.search(r'moderate',st): g='moderate'
    elif re.search(r'minimal',st): g='minimal'
    elif re.search(r'mild',st): g='mild'
    else: g='unspecified'
    fib=bool(re.search(r'fibros|cirrho|steatohepatitis|mallory|hyalin',t))
    rows.append(dict(donor=r.donor,grade=g,ord=ORD.get(g,np.nan),pct=p,fat=fat,neg=neg,fib_or_nash=fib,note=r.SMPTHNTS))
d=pd.DataFrame(rows)
d.to_csv(f'{L}/note_grades.tsv',sep='\t',index=False)
print(d.grade.value_counts().to_string())
print('graded (ord not NaN):',d.ord.notna().sum(),' stated pct:',d.pct.notna().sum(), ' pct>=33:',(d.pct>33).sum())
print('fib/NASH/cirrhosis/Mallory mention: all',d.fib_or_nash.sum(),'; among graded >= moderate',(d.fib_or_nash&(d.ord>=3)).sum())
print('UNSPEC',d[d.grade=='unspecified'].note.tolist())
for g in ['moderate','severe','minimal','none']:
    print(g, d[d.grade==g].note.str.slice(0,90).tolist()[:12])
