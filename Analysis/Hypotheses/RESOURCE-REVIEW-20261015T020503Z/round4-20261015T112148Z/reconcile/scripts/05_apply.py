"""Simulate applying the chosen edit set in memory (no repository file is written).
Decided variants: title clause removed (H52 untestable / H54 supports), H90 untestable,
E334 dropped, Methods R16-07 dropped where F1 R01-01 supersedes it.
Writes reconcile/applied/<path> copies and apply_log.json inside reconcile/ only."""
import json, os, re
ROOT = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas'
OUT = os.path.join(os.getcwd(), 'applied')
E = {e['id']: e for e in json.load(open('all_edits_located.json'))}
DROP = {
 'E188': 'Claims title kept-row variant', 'E190': 'Claims title kept-row variant', 'E192': 'Claims PAPER kept-row variant',
 'E323': 'F4abcd title variant not matching decided row', 'E324': 'same', 'E325': 'same', 'E326': 'same', 'E328': 'same',
 'E327': 'F4abcd title variant for decided row; Claims E187 chosen instead (also fixes clause 1)',
 'E067': 'H90 pass variant; H90 untestable', 'E070': 'H90 pass variant; H90 untestable',
 'E334': 'F2d placeholder note', 'E101': 'Methods R16-07 superseded by F1 R01-01', 'E103': 'same', 'E104': 'same',
}
# overlapping pairs: primary applied by span, secondary retried by string match after
PRIMARY_FIRST = {'E175': 'E019', 'E176': 'E020', 'E147': 'E237', 'E211': 'E282', 'E068': 'E293', 'E283': 'E078',
                 'E259': 'E084', 'E264': 'E085', 'E159': 'E065', 'E195': 'E213', 'E196': 'E214', 'E059': 'E321', 'E061': 'E153',
                 'E128': None}
log = []
byfile = {}
for i, e in E.items():
    if i in DROP:
        log.append(dict(id=i, action='dropped', why=DROP[i])); continue
    byfile.setdefault(e['file'], []).append(e)
applied_spans = {}
for f, es in byfile.items():
    p = os.path.join(ROOT, f)
    if not os.path.exists(p):
        for e in es: log.append(dict(id=e['id'], action='file_missing')); continue
    txt = open(p).read()
    first = [e for e in es if e['id'] not in PRIMARY_FIRST and e['found'] >= 1]
    later = [e for e in es if e['id'] in PRIMARY_FIRST or e['found'] < 1]
    # span edits, descending start; first occurrence only
    first.sort(key=lambda e: -e['spans'][0][0])
    last_start = len(txt) + 1
    marks = []
    for e in first:
        s, t = e['spans'][0][0], e['spans'][0][1]
        if t > last_start:
            later.append(e); continue
        txt = txt[:s] + e['new_text'] + txt[t:]
        last_start = s
        log.append(dict(id=e['id'], action='applied_span', file=f))
    for e in later:
        o = e['old_text']
        if o and o in txt:
            txt = txt.replace(o, e['new_text'], 1)
            log.append(dict(id=e['id'], action='applied_after_primary', file=f))
        else:
            log.append(dict(id=e['id'], action='NOT_APPLIED_old_text_gone', file=f, primary=PRIMARY_FIRST.get(e['id'])))
    q = os.path.join(OUT, f)
    os.makedirs(os.path.dirname(q), exist_ok=True)
    open(q, 'w').write(txt)
json.dump(log, open('apply_log.json', 'w'), indent=1)
from collections import Counter
print(Counter(l['action'] for l in log))
for l in log:
    if l['action'].startswith('NOT'): print(l)
