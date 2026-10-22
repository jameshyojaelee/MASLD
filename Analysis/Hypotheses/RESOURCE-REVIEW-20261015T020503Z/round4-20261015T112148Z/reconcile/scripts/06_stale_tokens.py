"""Numbers removed by an applied edit (in old_text, absent from new_text) that still occur
elsewhere in the post-application corpus. Reads reconcile/applied/ copies for edited files."""
import json, os, re, glob
ROOT = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas'
APP = os.path.join(os.getcwd(), 'applied')
E = {e['id']: e for e in json.load(open('all_edits_located.json'))}
log = json.load(open('apply_log.json'))
applied = [l['id'] for l in log if l['action'].startswith('applied')]
corpus = sorted(set(glob.glob(ROOT + '/docs/manuscript/**/*.md', recursive=True)) |
                {ROOT + '/docs/' + x for x in ['PAPER.md', 'RESULTS.md', 'STATUS.md', 'README.md', 'ROADMAP.md']} |
                set(glob.glob(ROOT + '/figures/main/*/README.md')))
corpus = [c for c in corpus if '/archive' not in c and os.path.exists(c)]
text = {}
for c in corpus:
    rel = os.path.relpath(c, ROOT)
    a = os.path.join(APP, rel)
    text[rel] = open(a if os.path.exists(a) else c).read()
NUM = re.compile(r'(?<![\w.\-−/])[-−]?\d[\d,]*(?:\.\d+)?%?(?![\w])')
def toks(s):
    out = set()
    for m in NUM.finditer(s):
        t = m.group(0).rstrip(',')
        d = re.sub(r'\D', '', t)
        if re.fullmatch(r'(19|20)\d\d', t): continue
        if len(d) >= 4 or ('.' in t and len(t.split('.')[-1].rstrip('%')) >= 3) or (',' in t):
            out.add(t.lstrip('-−'))
    return out
rows = []
for i in applied:
    e = E[i]
    gone = toks(e['old_text']) - toks(e['new_text'])
    for t in sorted(gone):
        pat = re.compile(r'(?<![\d.,])' + re.escape(t) + r'(?![\d])')
        for rel, tx in text.items():
            for m in pat.finditer(tx):
                ln = tx.count('\n', 0, m.start()) + 1
                line = tx.split('\n')[ln - 1]
                rows.append(dict(edit=i, ws=e['ws'], fid=e['fid'], edit_file=e['file'], token=t, hit_file=rel, hit_line=ln, context=line.strip()[:300]))
json.dump(rows, open('stale_hits.json', 'w'), indent=1)
from collections import Counter
print(len(rows), 'hits;', len({(r['token']) for r in rows}), 'tokens')
c = Counter((r['token'], r['ws'], r['fid']) for r in rows)
for k, v in sorted(c.items(), key=lambda x: -x[1]): print(v, k)
