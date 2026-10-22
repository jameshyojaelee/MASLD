import json, re
ROOT='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
RC=ROOT+'Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/'
mine={'R04-01','R04-02','R04-03','R04-04','R04-05','R04-07','R04-10','R10-01','R13-04','G1-01','G1-03'}
d=[e for e in json.load(open(RC+'repair/F4abcd/text_edits/structured_entries_repair.json')) if e['fid'] in mine]
# spot phrases copied from the task text
spot={'R04-01':['Stromal ECM is no longer among the rows shown in Figure\n4C','precision-weighted average, not one stage effect'],
 'R04-02':['3.6 for Ductular injury and 2.8 for Stromal ECM (variance ratios about 13 and\n8)'],
 'R04-04':['cholangiocyte Inflammatory / oxidative (SOD2) is first at 0.00307'],
 'R04-05':['In the T-cell fits, 17 of\n53 donors have fewer than 20 T cells (minimum 1)'],
 'R04-07':['carrying 0.64 of the frozen\nprogram weight'],
 'R10-01':['(empirical P=0.871, 200 draws)'],
 'G1-03':['smallest T-cell permutation p 0.1947, q 0.520']}
rows=[]
for e in d:
    alltxt=' '.join(t['new_text'] for t in e['text_edits'])
    for s in spot.get(e['fid'],[]): assert s in alltxt, (e['fid'],s)
    for t in e['text_edits']:
        txt=open(ROOT+t['file']).read()
        lines=txt.split('\n')
        n=txt.count(t['old_text'])
        a,b=(t['line'].split('-')+[t['line']])[:2]; a=int(a); b=int(b)
        seg='\n'.join(lines[a-1:b])
        at_line = t['old_text'] in seg or seg in t['old_text']
        new=t['new_text']
        em = '—' in new or ' -- ' in new
        tags = re.findall(r'\[(Inference|Unverified|Speculation|Estimate|Contains[^\]]*|Observed[^\]]*)\]', new)
        already = new in txt
        rows.append((e['fid'],t['file'],t['line'],n,at_line,already,em,tags))
        print(e['fid'],t['file'],t['line'],'old_count=',n,'at_line=',at_line,'new_present=',already,'emdash=',em,'tags=',tags)
