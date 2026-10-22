import json,sys
d=json.load(open(sys.argv[1]))
for a in d['articles']:
    ids=a.get('identifiers',{})
    j=a.get('journal'); j=j.get('title') if isinstance(j,dict) else j
    pd=a.get('publication_date') or {}
    print(ids.get('pmid') if isinstance(ids,dict) else ids, '|', j, '|', pd.get('year'), pd.get('month'), '|', a.get('doi'))
    print('  T:', a.get('title'))
    ab=a.get('abstract') or ''
    if isinstance(ab,dict): ab=' '.join(str(v) for v in ab.values())
    print('  A:', str(ab)[:int(sys.argv[2]) if len(sys.argv)>2 else 600].replace('\n',' '))
