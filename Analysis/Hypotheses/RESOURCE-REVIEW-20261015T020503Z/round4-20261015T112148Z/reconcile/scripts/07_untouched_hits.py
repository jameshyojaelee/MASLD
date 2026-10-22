"""Keep stale hits that sit in text no applied edit touched; report the ORIGINAL file:line.
Also grep a curated list of short stale tokens/phrases taken from the fix-agent cross-figure notes."""
import json, os, re, glob, sys
ROOT = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas'
E = {e['id']: e for e in json.load(open('all_edits_located.json'))}
log = json.load(open('apply_log.json'))
applied = [l['id'] for l in log if l['action'].startswith('applied')]
# touched original line ranges per file
touched = {}
for i in applied:
    e = E[i]
    if e['spans']:
        touched.setdefault(e['file'], []).append((e['spans'][0][2], e['spans'][0][3]))
corpus = sorted(set(glob.glob(ROOT + '/docs/manuscript/**/*.md', recursive=True)) |
                {ROOT + '/docs/' + x for x in ['PAPER.md', 'RESULTS.md', 'STATUS.md', 'README.md', 'ROADMAP.md']} |
                set(glob.glob(ROOT + '/figures/main/*/README.md')))
corpus = [c for c in corpus if '/reviews/' not in c and '/archive' not in c and os.path.exists(c)]
orig = {os.path.relpath(c, ROOT): open(c).read().split('\n') for c in corpus}
def untouched(rel, ln):
    return not any(a <= ln <= b for a, b in touched.get(rel, []))
CURATED = json.load(open(sys.argv[1])) if len(sys.argv) > 1 else []
rows = []
for item in CURATED:
    pat = re.compile(item['regex'])
    for rel, lines in orig.items():
        for n, line in enumerate(lines, 1):
            if pat.search(line):
                rows.append(dict(tag=item['tag'], file=rel, line=n, untouched=untouched(rel, n), text=line.strip()[:260]))
json.dump(rows, open('curated_hits.json', 'w'), indent=1)
for r in rows:
    if r['untouched']:
        print(f"{r['tag']:<28} {r['file']}:{r['line']}  {r['text'][:200]}")
