import sys, glob, json, xml.etree.ElementTree as ET
pm = {}
for f in glob.glob('sources/batch*.xml') + glob.glob('sources_v2/pubmed/*.xml'):
    try: root = ET.parse(f).getroot()
    except Exception: continue
    for a in root.iter('PubmedArticle'):
        p = a.findtext('.//PMID')
        ab = ' '.join(''.join(t.itertext()) for t in a.iter('AbstractText'))
        ids = {i.get('IdType'): i.text for i in a.iter('ArticleId')}
        pm[p] = (ab, ids.get('doi'))
for p in sys.argv[1:]:
    ab, doi = pm.get(p, ('MISSING', None))
    print('==', p, doi); print(ab); print()
