import json,re,sys
R="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
d=json.load(open(sys.argv[1]))
tags=re.compile(r"\[(Inference|Unverified|Speculation|Estimate|Contains Unverified Content)[^\]]*\]")
for fid,e in d.items():
  for t in (e if isinstance(e,list) else e.get("text_edits",[])):
    txt=open(R+t['file']).read()
    n=txt.count(t['old_text'])
    line=txt[:txt.find(t['old_text'])].count('\n')+1 if n else None
    new=t['new_text']
    probs=[]
    if '—' in new: probs.append('em dash')
    if '–' in new: probs.append('en dash')
    if tags.search(new): probs.append('evidence tag')
    for w in ['progression','causal','proves','demonstrat','ensures','novel']:
      if w in new.lower(): probs.append('word:'+w)
    print(f"{fid}\t{t['file']}\tclaimed {t['line']}\tfound_count={n}\tline={line}\t{probs}")
