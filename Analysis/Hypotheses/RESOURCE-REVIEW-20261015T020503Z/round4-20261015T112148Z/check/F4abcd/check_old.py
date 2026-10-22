import sys
sys.path.insert(0, sys.argv[1])
from edits_old import EDITS
ROOT='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
for fid,f,line,old in EDITS:
    txt=open(ROOT+f).read()
    lines=txt.split('\n')
    a,b=(line.split('-')+[line])[:2]
    a,b=int(a),int(b)
    seg='\n'.join(lines[a-1:b])
    n=txt.count(old)
    at = old in seg or seg in old and seg==old
    # find line of first occurrence
    pos=txt.find(old); ln = txt[:pos].count('\n')+1 if pos>=0 else None
    endln = ln+old.count('\n') if ln else None
    print(f"{fid}\t{f}\t{line}\tcount={n}\tfound_at={ln}-{endln}\tline_match={'yes' if ln==a and endln==b else 'NO'}")
