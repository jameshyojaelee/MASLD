import json, os
ROOT = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas'
def load(f):
    d = json.load(open(f))
    return d['final_edits'] if 'final_edits' in d else d['edits']
gen = load('staged/edits_gen.json')
others = []
for w in ('bulk', 'cat', 'sc'):
    for e in load(f'staged/final_{w}.json'):
        others.append((w, e))
src_cache = {}
def span(e):
    s = src_cache.setdefault(e['file'], open(os.path.join(ROOT, e['file'])).read())
    i = s.find(e['old']); return i, i + len(e['old'])
for g in gen:
    a = span(g)
    for w, o in others:
        if o['file'] != g['file']: continue
        b = span(o)
        if b[0] < a[1] and a[0] < b[1]:
            print('OVERLAP', g['id'], a, w, o['id'], b)
            print('  GEN OLD:', g['old'][:300]); print('  GEN NEW:', g['new'][:600])
            print('  OTH OLD:', o['old'][:300]); print('  OTH NEW:', o['new'][:600])
for i in ('cat-12','cat-78','cat-07','cat-08'):
    print(i, 'in final_cat:', any(o['id']==i for w,o in others), 'dropped:', [x for x in json.load(open('staged/final_cat.json')).get('dropped',[]) if x['id']==i])
