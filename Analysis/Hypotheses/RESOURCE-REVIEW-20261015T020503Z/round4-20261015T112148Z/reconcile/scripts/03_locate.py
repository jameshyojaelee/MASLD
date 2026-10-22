"""Locate every proposed old_text in the current file; record spans and line ranges."""
import json, os, re
ROOT = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas'
R4 = os.path.dirname(os.getcwd())
A = json.load(open('agent_outputs.json'))
edits = []
for k, v in A.items():
    if not k.startswith('fix:') or not v:
        continue
    ws = k[4:]
    for e in v.get('entries', []):
        for i, t in enumerate(e.get('text_edits') or []):
            edits.append(dict(ws=ws, fid=e['fid'], status=e.get('status'), file=t.get('file'), line=t.get('line'),
                              claim_level=t.get('claim_level'), old_text=t.get('old_text') or '', new_text=t.get('new_text') or '', src='structured_output'))
# F2d: structured output abbreviates some edits; replace with the saved full-text file
f2d = json.load(open(R4 + '/ws/F2d/text_edits.json'))
def abbreviated(t):
    return 'text_edits.json' in t['new_text'] or 'text_edits.json' in t['old_text']
full_used = set()
for e in edits:
    if e['ws'] == 'F2d' and abbreviated(e):
        cands = [i for i, f in enumerate(f2d) if f['fid'] == e['fid'] and f['file'] == e['file'] and (e['old_text'][:40] in f['old_text'] or str(f.get('line')) == str(e['line']))]
        if cands:
            f = f2d[cands[0]]; full_used.add(cands[0])
            e.update(old_text=f['old_text'], new_text=f['new_text'], src='ws/F2d/text_edits.json')
        else:
            e['src'] = 'ABBREVIATED_NO_MATCH'
# F2d file edits not in structured output
so_keys = {(e['file'], e['old_text']) for e in edits if e['ws'] == 'F2d'}
for i, f in enumerate(f2d):
    if (f['file'], f['old_text']) not in so_keys:
        edits.append(dict(ws='F2d', fid=f['fid'], status='file_only', file=f['file'], line=f.get('line'), claim_level=f.get('claim_level'),
                          old_text=f['old_text'], new_text=f['new_text'], src='ws/F2d/text_edits.json (not in structured output)'))
cache = {}
for n, e in enumerate(edits):
    e['id'] = 'E%03d' % n
    p = os.path.join(ROOT, e['file'])
    if p not in cache:
        cache[p] = open(p).read() if os.path.exists(p) else None
    txt = cache[p]
    e['found'] = 0; e['spans'] = []
    if txt is None:
        e['found'] = -1; continue
    o = e['old_text']
    if not o:
        continue
    start = 0
    while True:
        i = txt.find(o, start)
        if i < 0: break
        l0 = txt.count('\n', 0, i) + 1
        l1 = txt.count('\n', 0, i + len(o)) + 1
        e['spans'].append([i, i + len(o), l0, l1])
        start = i + 1
    e['found'] = len(e['spans'])
    if e['found'] == 0:
        # whitespace-normalised search for diagnostics
        norm = lambda s: re.sub(r'\s+', ' ', s)
        e['found_ws_normalised'] = norm(o) in norm(txt)
json.dump(edits, open('all_edits_located.json', 'w'), indent=1)
from collections import Counter
print(len(edits), Counter(e['found'] for e in edits))
for e in edits:
    if e['found'] != 1:
        print(e['id'], e['ws'], e['fid'], e['file'], e['line'], 'found=', e['found'], 'wsnorm=', e.get('found_ws_normalised'), repr(e['old_text'][:90]))
