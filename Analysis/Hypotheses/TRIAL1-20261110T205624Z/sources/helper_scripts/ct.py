import json,sys,urllib.request,time,os
for n in sys.argv[1].split(','):
    p='raw/%s.json'%n
    if not os.path.exists(p):
        open(p,'wb').write(urllib.request.urlopen('https://clinicaltrials.gov/api/v2/studies/'+n).read()); time.sleep(0.3)
    d=json.load(open(p))['protocolSection']
    s=d['statusModule']; de=d['designModule']; c=d.get('conditionsModule',{}).get('conditions',[])
    print(n,'|',d['identificationModule'].get('acronym',''),'|',de.get('phases'),'| n=',de.get('enrollmentInfo',{}).get('count'),de.get('enrollmentInfo',{}).get('type'),'|',s.get('overallStatus'),'|',d['identificationModule']['briefTitle'][:90],'|',c[:3], '| results' if json.load(open(p)).get('hasResults') else '')
