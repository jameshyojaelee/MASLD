"""Panel-deep install (user approved decisions 1-9 on 2026-10-17: 'yes to all'). Archives every overwritten or retired file first."""
import csv, shutil, os, sys, hashlib, json, re
ROOT='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'; os.chdir(ROOT)
P='figures/candidates/panel-deep-fixes-20261017T035633Z/'; INST=sys.argv[1]; ARC=INST+'/archive/'
log=[]
def sha(p): return hashlib.sha256(open(p,'rb').read()).hexdigest()
def archive(t):
    if os.path.exists(t):
        d=ARC+t; assert not os.path.exists(d),d; os.makedirs(os.path.dirname(d),exist_ok=True); shutil.copy2(t,d)
def put(src,t,kind):
    archive(t); os.makedirs(os.path.dirname(t),exist_ok=True); shutil.copy2(src,t)
    log.append({'kind':kind,'source':src,'target':t,'sha256':sha(t)})
for r in csv.DictReader(open(P+'INSTALL_PLAN.tsv'),delimiter='\t'):
    if r['status']!='redrawn': continue
    assert os.path.exists(r['candidate_pdf']) and os.path.exists(r['install_target']),r['panel']
    put(r['candidate_pdf'],r['install_target'],'panel:'+r['panel']+':decision'+(r['decision'] or '-'))
    st=r['source_table']; srcs=[st.split(' (+')[0]]+[os.path.join(os.path.dirname(st.split(' (+')[0]),x) for x in re.findall(r'\(\+ ([^)]+)\)',st)]
    dd=os.path.join(os.path.dirname(r['install_target']),'data','panel_deep_20261017')
    for s in srcs:
        assert os.path.exists(s),s; put(s,os.path.join(dd,os.path.basename(s)),'source_table:'+r['panel'])
# retire the stale August 5E selection table (13 rows; reads neither protein table)
old='figures/main/fig5_molecular_context/panels/data/fig5e_multimodal_program_summary_selection.tsv'
if os.path.exists(old):
    archive(old); d='figures/main/fig5_molecular_context/panels/_archive/panel-deep-'+os.path.basename(INST)+'/'
    os.makedirs(d,exist_ok=True); shutil.move(old,d+os.path.basename(old)); log.append({'kind':'retire','target':old,'moved_to':d})
json.dump(log,open(INST+'/install_log.json','w'),indent=1)
print(len(log),'actions;',sum(x['kind'].startswith('panel') for x in log),'panels;',sum(x['kind'].startswith('source') for x in log),'source tables')
