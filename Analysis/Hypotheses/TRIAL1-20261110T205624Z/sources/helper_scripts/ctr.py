import json,sys,urllib.request,os
n=sys.argv[1]; p='raw/%s.json'%n
if not os.path.exists(p): open(p,'wb').write(urllib.request.urlopen('https://clinicaltrials.gov/api/v2/studies/'+n).read())
d=json.load(open(p)); ps=d['protocolSection']
print(n, ps['designModule'].get('designInfo',{}).get('allocation'), ps['designModule'].get('designInfo',{}).get('maskingInfo',{}).get('masking'))
kw=sys.argv[2].lower().split(',') if len(sys.argv)>2 else None
rs=d.get('resultsSection',{})
for om in rs.get('outcomeMeasuresModule',{}).get('outcomeMeasures',[]):
    t=om['title']
    if kw and not any(k in t.lower() for k in kw): continue
    print('\n#',om.get('type'),'|',t,'|',om.get('timeFrame'),'|',om.get('unitOfMeasure'),'|',om.get('paramType'))
    g={x['id']:x['title'] for x in om.get('groups',[])}
    for c in om.get('classes',[]):
        for cat in c.get('categories',[]):
            for m in cat.get('measurements',[]):
                print('  ',c.get('title',''),cat.get('title',''),g.get(m['groupId']),m.get('value'),m.get('spread',''),m.get('lowerLimit',''),m.get('upperLimit',''))
    for a in om.get('analyses',[]):
        print('  AN:',[g.get(x) for x in a.get('groupIds',[])],a.get('pValue'),a.get('statisticalMethod'),a.get('paramType'),a.get('paramValue'),a.get('ciLowerLimit'),a.get('ciUpperLimit'),a.get('nonInferiorityType'),(a.get('pValueComment') or '')[:150])
