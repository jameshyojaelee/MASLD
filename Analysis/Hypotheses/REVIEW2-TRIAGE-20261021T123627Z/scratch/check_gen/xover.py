import json, glob, os
ROOT = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas'
spans = {}
for f in sorted(glob.glob('staged/edits_*.json') + glob.glob('staged/final_*.json') + glob.glob('staged/main_*.json')):
    d = json.load(open(f))
    eds = d['edits'] if isinstance(d, dict) and 'edits' in d else d.get('final_edits', d) if isinstance(d, dict) else d
    for e in eds:
        src = open(os.path.join(ROOT, e['file'])).read()
        i = src.find(e['old'])
        spans.setdefault(e['file'], []).append((i, i + len(e['old']), os.path.basename(f), e['id'], src.count(e['old'])))
for fn, ss in spans.items():
    ss.sort()
    for a, b in zip(ss, ss[1:]):
        if b[0] < a[1] and a[2] != b[2]:
            print('CROSS OVERLAP', fn, a, b)
        elif b[0] < a[1]:
            print('SAME-WRITER OVERLAP', fn, a, b)
    for s in ss:
        if s[4] != 1: print('COUNT', fn, s)
print('done')
