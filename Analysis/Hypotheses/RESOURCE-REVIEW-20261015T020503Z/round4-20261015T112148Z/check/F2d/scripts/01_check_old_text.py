import json, os, re
ROOT='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas'
E=json.load(open(ROOT+'/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/ws/F2d/text_edits.json'))
out=[]
for i,e in enumerate(E):
    f=os.path.join(ROOT,e['file']); txt=open(f).read()
    old=e.get('old_text'); new=e.get('new_text','')
    if old is None:
        print(i,e['fid'],e['file'],'NO old_text', e.get('note')); continue
    n=txt.count(old)
    if n:
        pos=txt.index(old); ln=txt[:pos].count('\n')+1; lnend=ln+old.count('\n')
    else:
        ln=lnend=None
    em = ('—' in new)
    endash = '–' in new
    print(f"{i:2d} {e['fid']} {e['file']} claimed_line={e.get('line')} found={n} at {ln}-{lnend} emdash_new={em} endash_new={endash}")
    out.append(dict(i=i,fid=e['fid'],file=e['file'],claimed=e.get('line'),count=n,start=ln,end=lnend))
json.dump(out,open('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/F2d/out/old_text_check.json','w'),indent=1)
