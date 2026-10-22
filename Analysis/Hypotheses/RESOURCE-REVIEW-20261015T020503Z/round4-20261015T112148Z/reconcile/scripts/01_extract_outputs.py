"""Extract final StructuredOutput payloads of fix:* and check:* agents of the round-4 workflow."""
import json, glob, os
D = '/gpfs/commons/home/jameslee/.claude/projects/-gpfs-commons-groups-sanjana-lab-Cas13-MASLD-library-design/39a23967-7721-4ff8-86d5-d2d5894e47f5/subagents/workflows/wf_c93b4f76-46f'
OUT = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/reconcile'
res = {}
for m in sorted(glob.glob(D + '/*.meta.json')):
    desc = json.load(open(m))['description']
    if not (desc.startswith('fix:') or desc.startswith('check:')):
        continue
    j = m.replace('.meta.json', '.jsonl')
    last = None
    for line in open(j):
        try:
            rec = json.loads(line)
        except Exception:
            continue
        msg = rec.get('message') or {}
        content = msg.get('content') if isinstance(msg, dict) else None
        if isinstance(content, list):
            for c in content:
                if isinstance(c, dict) and c.get('type') == 'tool_use' and c.get('name') == 'StructuredOutput':
                    last = c.get('input')
    key = desc
    if key in res:
        key = desc + '#' + os.path.basename(j)[:15]
    res[key] = last
    print(desc, 'found' if last else 'NONE', len(json.dumps(last)) if last else 0)
json.dump(res, open(OUT + '/agent_outputs.json', 'w'), indent=1)
