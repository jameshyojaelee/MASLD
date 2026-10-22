"""Check F4efgh repair text edits against the current files: verbatim, unique, no overlap, no em dash, no evidence tags."""
import glob, json, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import new_edits

ROOT = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
R4 = ROOT + 'Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/'
r1 = json.load(open(R4 + 'repair_inputs/F4efgh.json'))
edits = list(new_edits.E)
# carry unchanged first-round edits for R05-07 and R05-09
for e in r1['entries']:
    if e['fid'] in ('R05-07', 'R05-09'):
        for t in e.get('text_edits') or []:
            edits.append(dict(fid=e['fid'], **{k: t.get(k) for k in ('file', 'line', 'old_text', 'new_text', 'claim_level')}))

bad = 0
spans = {}
for i, t in enumerate(edits):
    txt = open(ROOT + t['file']).read()
    n = txt.count(t['old_text'])
    if n != 1:
        bad += 1
        print('FAIL count=%d' % n, t['fid'], t['file'], t['line'], repr(t['old_text'][:70]))
        continue
    s = txt.index(t['old_text'])
    ln = txt[:s].count('\n') + 1
    spans.setdefault(t['file'], []).append((s, s + len(t['old_text']), t['fid'], ln, i))
    if '—' in t['new_text']:
        bad += 1; print('EM DASH', t['fid'], t['file'])
    if re.search(r'\[(Inference|Unverified|Speculation|Estimate|Contains)', t['new_text']):
        bad += 1; print('TAG', t['fid'], t['file'])
    if str(ln) != str(t['line']).split('-')[0]:
        print('line note: %s %s stated %s actual %d' % (t['fid'], t['file'], t['line'], ln))
for f, sp in spans.items():
    sp.sort()
    for a, b in zip(sp, sp[1:]):
        if b[0] < a[1]:
            bad += 1; print('OVERLAP (mine)', f, a[2], a[3], b[2], b[3])
# other workstreams' round-1 edits
others = []
for p in sorted(glob.glob(R4 + 'repair_inputs/*.json')):
    w = os.path.basename(p)[:-5]
    if w in ('F4efgh', 'reconcile_conflicts'):
        continue
    d = json.load(open(p))
    for e in d.get('entries', []):
        for t in e.get('text_edits') or []:
            if t['file'] in spans:
                txt = open(ROOT + t['file']).read()
                if txt.count(t['old_text']) != 1:
                    continue
                s = txt.index(t['old_text'])
                for (a0, a1, fid, ln, i) in spans[t['file']]:
                    if s < a1 and a0 < s + len(t['old_text']):
                        print('OVERLAP other', w, e['fid'], t['file'], 'line', txt[:s].count('\n') + 1, 'with', fid, ln)
                    elif abs(s - a1) < 3 or abs(a0 - (s + len(t['old_text']))) < 3:
                        print('ADJACENT other', w, e['fid'], t['file'], 'with', fid, ln)
print('edits', len(edits), 'bad', bad)
json.dump(edits, open(os.path.dirname(os.path.abspath(__file__)) + '/tmp/all_edits.json', 'w'), indent=1)
