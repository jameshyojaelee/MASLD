import json
A = json.load(open('agent_outputs.json'))
edits = []
cross = []
for k, v in A.items():
    if not k.startswith('fix:') or not v:
        continue
    ws = v.get('ws') or k[4:]
    for e in v.get('entries', []):
        for t in e.get('text_edits') or []:
            edits.append(dict(ws=ws, fid=e['fid'], status=e.get('status'), file=t.get('file'), line=t.get('line'), claim_level=t.get('claim_level'), old_text=t.get('old_text'), new_text=t.get('new_text')))
    for c in v.get('cross_figure_numbers') or []:
        cross.append(dict(ws=ws, text=c))
json.dump(edits, open('all_edits.json', 'w'), indent=1)
json.dump(cross, open('cross_numbers.json', 'w'), indent=1)
print(len(edits), 'edits', len(cross), 'cross')
from collections import Counter
print(Counter(e['ws'] for e in edits))
print(Counter(e['file'] for e in edits).most_common())
# statuses
for k, v in A.items():
    if k.startswith('fix:') and v:
        print(k, Counter(e['status'] for e in v['entries']))
