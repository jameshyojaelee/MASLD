import json, re, os
root='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas'
E=json.load(open('edits_input.json'))
for f in E:
    for i,e in enumerate(f['text_edits']):
        p=os.path.join(root,e['file']); t=open(p).read()
        n_old=t.count(e['old_text']); n_new=t.count(e['new_text'])
        # line location
        idx=t.find(e['old_text']); ln = t[:idx].count('\n')+1 if idx>=0 else None
        nt=e['new_text']
        issues=[]
        if '—' in nt: issues.append('em dash')
        if re.search(r'\[(Inference|Unverified|Speculation|Estimate|Contains)', nt): issues.append('evidence tag')
        print(f"{f['fid']} #{i} {e['file']} claimed {e['line']} found_line={ln} old_count={n_old} new_count={n_new} {issues}")
