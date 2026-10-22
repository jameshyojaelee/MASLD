# Re-apply the 7A producer's counting rule (00_extract_panel_data.py:591-616) to the current audit files
import json, glob, collections, os
base='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/config/artifacts/models'
states=collections.Counter()
for p in sorted(glob.glob(base+'/*/exposure_audit.json')):
    aj=json.load(open(p)); found=[]
    for block in ("checkpoint_findings","checkpoint_finding"):
        v=aj.get(block)
        if isinstance(v,dict) and "exposure_state" in v: found.append(v["exposure_state"])
        elif isinstance(v,dict):
            for _k,_v in v.items():
                if isinstance(_v,dict) and _v.get("exposure_state"): found.append(_v["exposure_state"])
    if not found: states["no finding"]+=1
    else:
        for s in found:
            if s=="unknown": states["unknown"]+=1
            elif s.startswith("target_label_unexposed"): states["unexposed"]+=1
            elif "seen" in s: states["seen"]+=1
            else: states["other"]+=1
print(dict(states), 'total', sum(states.values()), 'findings', sum(v for k,v in states.items() if k!='no finding'))
