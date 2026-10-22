"""Step 2: 7B/7D/7E producer-value panels, then the 137 merged text edits. Archives first (archive/step2/)."""
import json, shutil, os, sys, hashlib
ROOT='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'; os.chdir(ROOT)
INST=sys.argv[1]; TM=sys.argv[2]; ARC=INST+'/archive/step2/'
C='figures/candidates/panel-deep-fixes-20261017T035633Z/fig7/'; T='figures/main/fig7_model_framework/'
log=[]
def sha(p): return hashlib.sha256(open(p,'rb').read()).hexdigest()
def archive(t):
    d=ARC+t; assert not os.path.exists(d),d; os.makedirs(os.path.dirname(d),exist_ok=True); shutil.copy2(t,d)
def put(src,t,kind):
    if os.path.exists(t): archive(t)
    os.makedirs(os.path.dirname(t),exist_ok=True); shutil.copy2(src,t); log.append({'kind':kind,'source':src,'target':t,'sha256':sha(t)})
for p,stem in [('7B','fig7B_gain_over_task_native_baseline'),('7D','fig7D_measured_bounds_and_join_structure'),('7E','fig7E_verdict_vocabulary_and_refusal_thresholds')]:
    put(C+'panels/'+stem+'.producer.pdf', T+stem+'.pdf', 'panel_producer:'+p)
    put(C+'source_tables/fig'+p+'.producer.tsv', T+'data/panel_deep_20261017/fig'+p+'.tsv', 'source_table_producer:'+p)
E=json.load(open(TM+'/edits.json')); files={}
for e in E: files.setdefault(e['file'],[]).append(e)
texts={}
for f,es in files.items():
    s=open(f).read()
    for e in es:
        assert s.count(e['old'])==1,(e['id'],f,s.count(e['old'])); s=s.replace(e['old'],e['new'])
    texts[f]=s
for f,s in texts.items():
    archive(f); open(f,'w').write(s); log.append({'kind':'text','target':f,'n_edits':len(files[f]),'sha256':sha(f)})
json.dump(log,open(INST+'/install_step2_log.json','w'),indent=1)
print(len(log),'actions;',sum(len(v) for v in files.values()),'edits in',len(files),'files')
