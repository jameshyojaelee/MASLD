import json, glob, collections
base='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/config/artifacts/models'
files=sorted(glob.glob(base+'/*/exposure_audit.json'))
states=collections.Counter(); nf=0; none=0
for f in files:
    d=json.load(open(f)); found=[]
    cf=d.get('checkpoint_findings')
    if isinstance(cf,dict): found+= [v for v in cf.values()]
    c1=d.get('checkpoint_finding')
    if isinstance(c1,dict):
        found.append(c1)
    if found:
        for v in found:
            nf+=1; states[v.get('exposure_state') if isinstance(v,dict) else 'non-dict']+=1
    else: none+=1
print('files',len(files),'findings',nf,'files without any finding',none,'entries',nf+none)
for k,v in states.most_common(): print(' ',v,k)
