import json, os
ROOT = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
edits = json.load(open('tmp/all_edits.json'))
files = {}
for t in edits:
    files.setdefault(t['file'], open(ROOT + t['file']).read())
out = []
for t in sorted(edits, key=lambda t: (t['file'], files[t['file']].find(t['old_text']))):
    txt = files[t['file']]
    s = txt.index(t['old_text'])
    a = txt.rfind('\n', 0, max(0, txt.rfind('\n', 0, s)))  # one line before
    e = s + len(t['old_text'])
    b = txt.find('\n', txt.find('\n', e) + 1)
    out.append('##### %s %s\n%s[[%s]]%s\n' % (t['fid'], t['file'], txt[a+1:s], t['new_text'], txt[e:b]))
open('tmp/preview.txt', 'w').write('\n'.join(out))
