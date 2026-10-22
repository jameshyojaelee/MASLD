import json, os
ROOT = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
OUT = 'tmp/applied/'
edits = json.load(open('tmp/all_edits.json'))
files = {}
for t in edits:
    files.setdefault(t['file'], open(ROOT + t['file']).read())
for t in edits:
    assert files[t['file']].count(t['old_text']) == 1, (t['fid'], t['file'])
    files[t['file']] = files[t['file']].replace(t['old_text'], t['new_text'])
for f, txt in files.items():
    p = OUT + f
    os.makedirs(os.path.dirname(p), exist_ok=True)
    open(p, 'w').write(txt)
print('applied', len(edits), 'to', len(files), 'files')
