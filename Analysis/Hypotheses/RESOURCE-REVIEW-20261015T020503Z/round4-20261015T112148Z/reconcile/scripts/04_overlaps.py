import json, itertools
E = json.load(open('all_edits_located.json'))
by = {}
for e in E:
    by.setdefault(e['file'], []).append(e)
rows = []
for f, es in by.items():
    for a, b in itertools.combinations(es, 2):
        for sa in a['spans']:
            for sb in b['spans']:
                if sa[0] < sb[1] and sb[0] < sa[1]:
                    same = (a['old_text'] == b['old_text'])
                    samenew = (a['new_text'] == b['new_text'])
                    rows.append(dict(file=f, a=a['id'], b=b['id'], a_ws=a['ws'], b_ws=b['ws'], a_fid=a['fid'], b_fid=b['fid'],
                                     a_line=a['line'], b_line=b['line'], a_lines=sa[2:], b_lines=sb[2:], same_old=same, same_new=samenew))
json.dump(rows, open('overlaps.json', 'w'), indent=1)
for r in rows:
    print(r['file'], r['a'], r['a_ws'], r['a_fid'], r['a_lines'], '|', r['b'], r['b_ws'], r['b_fid'], r['b_lines'], 'same_old' if r['same_old'] else '', 'same_new' if r['same_new'] else '')
print(len(rows))
