"""Fetch a ClinicalTrials.gov v2 record (fresh, saved to raw/) and print outcome results and analyses.
Usage: python3 ctr2.py NCTID [comma-separated title keywords] [--full]"""
import json, sys, urllib.request, os
D = os.path.dirname(os.path.abspath(__file__))
n = sys.argv[1]; p = os.path.join(D, 'raw', n + '.json')
if not os.path.exists(p):
    open(p, 'wb').write(urllib.request.urlopen('https://clinicaltrials.gov/api/v2/studies/' + n).read())
d = json.load(open(p)); ps = d['protocolSection']
print(n, ps['identificationModule'].get('acronym', ''), ps['designModule'].get('phases'),
      'enrol', ps['designModule'].get('enrollmentInfo'), 'hasResults', d.get('hasResults'))
kw = [k for k in (sys.argv[2].lower().split(',') if len(sys.argv) > 2 and sys.argv[2] != '-' else [])]
rs = d.get('resultsSection', {})
for om in rs.get('outcomeMeasuresModule', {}).get('outcomeMeasures', []):
    t = om['title']
    if kw and not any(k in t.lower() for k in kw): continue
    print('\n#', om.get('type'), '|', t, '|', om.get('timeFrame'), '|', om.get('unitOfMeasure'), '|', om.get('paramType'))
    if '--full' in sys.argv: print('  DESC:', om.get('description', '')[:1500]); print('  POP:', om.get('populationDescription', '')[:600])
    g = {x['id']: x['title'] for x in om.get('groups', [])}
    for den in om.get('denoms', []):
        print('  N:', [(g.get(c['groupId']), c['value']) for c in den.get('counts', [])])
    for c in om.get('classes', []):
        for cat in c.get('categories', []):
            for m in cat.get('measurements', []):
                print('  ', c.get('title', ''), cat.get('title', ''), '|', g.get(m['groupId']), m.get('value'), m.get('spread', ''), m.get('lowerLimit', ''), m.get('upperLimit', ''))
    for a in om.get('analyses', []):
        print('  AN:', [g.get(x) for x in a.get('groupIds', [])], 'p=', a.get('pValue'), a.get('statisticalMethod'), a.get('paramType'), a.get('paramValue'), a.get('ciPctValue'), a.get('ciLowerLimit'), a.get('ciUpperLimit'))
        for k in ('groupDescription', 'pValueComment', 'statisticalComment', 'estimateComment'):
            if a.get(k): print('     ', k, ':', a[k][:700])
