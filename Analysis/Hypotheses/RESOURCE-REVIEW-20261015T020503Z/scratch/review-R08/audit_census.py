import json, glob, os, collections
base='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/config/artifacts/models'
dirs=[d for d in sorted(os.listdir(base)) if os.path.isdir(os.path.join(base,d))]
files=sorted(glob.glob(base+'/*/exposure_audit.json'))
states=collections.Counter(); nfind=0; nofind=0; per=[]
for f in files:
    d=json.load(open(f)); cf=d.get('checkpoint_findings')
    if isinstance(cf,dict) and len(cf)>0:
        for k,v in cf.items():
            nfind+=1; states[(v.get('exposure_state') if isinstance(v,dict) else 'non-dict')]+=1
    else:
        nofind+=1
print('model dirs',len(dirs),'audit files',len(files))
print('checkpoint findings',nfind,'files with no finding',nofind,'total entries',nfind+nofind)
for k,v in states.most_common(): print(' ',k,v)
