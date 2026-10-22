import pandas as pd, re, numpy as np
L='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L02-01'
m=pd.read_csv(f'{L}/liver_rna_meta.tsv',sep='\t')
def pct(t):
    # take numbers immediately tied to % (ranges -> mean); '>80%','80+%' -> value
    t=t.replace('!0%','10%').replace('%60','60%')
    vals=[]
    for mm in re.finditer(r'(\d+)\s*(?:-|to)\s*(\d+)\s*\+?\s*%|([<>~]?)\s*(\d+)\s*\+?\s*%',t):
        if mm.group(1): vals.append((float(mm.group(1))+float(mm.group(2)))/2)
        else:
            v=float(mm.group(4)); 
            if mm.group(3)=='<': v=v/2
            vals.append(v)
    return vals
rows=[]
for _,r in m.iterrows():
    t=str(r.SMPTHNTS).lower()
    fat=bool(re.search(r'steato|fat',t))
    neg=bool(re.search(r'no significant steatosis|no steatosis|without steatosis|no major|no fat|no significant fat|no obvious steatosis',t))
    p=None
    if fat:
        # percentages attached to fat context: only if note mentions fat; exclude fibrous tissue %
        seg=re.split(r'[;.]',t)
        pv=[]
        for s in seg:
            if re.search(r'steato|fat',s): pv+=pct(s)
        if pv: p=max(pv)
    # grade
    if not fat: g='no_mention'
    elif neg: g='none'
    elif p is not None:
        g='minimal' if p<5 else ('mild' if p<=33 else ('moderate' if p<=66 else 'severe'))
        if p==5: g='mild'
    elif re.search(r'severe|marked|extensive|panlobular|>\s*80|90',t): g='severe'
    elif re.search(r'moderate',t) and re.search(r'mild[- ]to[- ]moderate|mild-moderate',t): g='mild_moderate'
    elif re.search(r'moderate',t) and re.search(r'steato|fat',t): g='moderate'
    elif re.search(r'minimal',t): g='minimal'
    elif re.search(r'mild',t): g='mild'
    elif re.search(r'diffuse',t): g='diffuse_unspecified'
    else: g='unspecified'
    fib=bool(re.search(r'fibros|cirrho|steatohepatitis|mallory|hyalin',t))
    rows.append(dict(donor=r.donor,grade=g,pct=p,fat=fat,neg=neg,fib_or_nash=fib,note=r.SMPTHNTS))
d=pd.DataFrame(rows)
d.to_csv(f'{L}/note_grades.tsv',sep='\t',index=False)
print(d.grade.value_counts().to_string())
print('stated pct n',d.pct.notna().sum())
print('fibrosis/NASH/cirrhosis/Mallory mention among fat-mention:',(d.fib_or_nash&d.fat&~d.neg).sum(),' among all:',d.fib_or_nash.sum())
print(d[d.grade=='unspecified'].note.tolist())
print(d[d.grade=='moderate'].note.head(30).tolist())
