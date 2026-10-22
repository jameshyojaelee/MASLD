"""Round-4 install (user approved A-S except O on 2026-10-16). Archives every overwritten file first."""
import json, csv, shutil, os, sys, hashlib
ROOT='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'; os.chdir(ROOT)
R4='Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/'
ST=R4+'staged-20261015T184532Z/'; INST=sys.argv[1]; ARC=INST+'/archive/'
log=[]
def sha(p): return hashlib.sha256(open(p,'rb').read()).hexdigest()
def archive(t):
    if os.path.exists(t):
        d=ARC+t; os.makedirs(os.path.dirname(d),exist_ok=True); shutil.copy2(t,d)
def put(src,t,kind):
    archive(t); os.makedirs(os.path.dirname(t),exist_ok=True); shutil.copy2(src,t)
    log.append({'kind':kind,'source':src,'target':t,'sha256':sha(t)})
# 1 conditional text fixes (C, F, G approved; fig4 note)
for e in json.load(open(R4+'wave4/conditional_text_fixes.json')):
    p=ST+'final/'+e['file']; s=open(p).read(); assert s.count(e['old_text'])==1,e['file']
    open(p,'w').write(s.replace(e['old_text'],e['new_text'])); log.append({'kind':'conditional_fix','target':e['file']})
# 2 text
for dp,_,fs in os.walk(ST+'final'):
    for f in fs:
        src=os.path.join(dp,f); t=os.path.relpath(src,ST+'final'); put(src,t,'text')
# 3 panels and sidecars
for r in csv.DictReader(open(R4+'wave4/INSTALL_PLAN.tsv'),delimiter='\t'):
    if r['item'] in ('A','O','Q','ALL'): continue
    src=r['source'] if os.path.exists(r['source']) else R4+r['source']
    assert os.path.exists(src),src
    put(src,r['target'],'panel:'+r['item'])
# 4 retire superseded Fig 1B summary (OR 0.757)
old='figures/main/fig1_atlas_overview/fig1b_interface_summary.tsv'
if os.path.exists(old):
    d='figures/main/fig1_atlas_overview/_archive/round4-'+os.path.basename(INST)+'/'; os.makedirs(d,exist_ok=True)
    archive(old); shutil.move(old,d+os.path.basename(old)); log.append({'kind':'retire','target':old,'moved_to':d})
# 5 item A validator markers
title=' '.join(open(ST+'final/docs/manuscript/draft/title.md').read().split('**')[1].split())
reps={'config/resource_paper_scope.config.json':[('A multimodal human MASLD Resource resolves the architecture of inherited risk and the coordination of multicellular remodeling',title)],
 'scripts/manuscript/validate_resource_scope.py':[('"The paper connects **inherited genetic variation** to **multicellular disease remodeling** in human MASLD."','"The paper records **liver-trait genetic evidence** and **disease-state remodeling** in human MASLD side by side for each gene."'),
   ('"| Catalog strict `tested_negative` calls | 0 |"','"| Catalog `tested_negative` calls | not assessed |"')]}
for p,rs in reps.items():
    archive(p); s=open(p).read()
    for a,b in rs: assert s.count(a)==1,(p,a); s=s.replace(a,b)
    open(p,'w').write(s); log.append({'kind':'validator_marker','target':p,'title':title})
json.dump(log,open(INST+'/install_log.json','w'),indent=1)
print(len(log),'actions;',sum(1 for x in log if x['kind']=='text'),'text,',sum(1 for x in log if x['kind'].startswith('panel')),'panels')
