import json,sys,urllib.request,urllib.parse
q=sys.argv[1]
u='https://clinicaltrials.gov/api/v2/studies?pageSize=30&fields=NCTId,BriefTitle,Phase,OverallStatus,EnrollmentCount,Condition,HasResults&'+q
d=json.load(urllib.request.urlopen(u))
for s in d['studies']:
    p=s['protocolSection']
    print(p['identificationModule']['nctId'],p.get('designModule',{}).get('phases'),p.get('designModule',{}).get('enrollmentInfo',{}).get('count'),p['statusModule']['overallStatus'],'|',p['identificationModule']['briefTitle'][:110], '| R' if s.get('hasResults') else '')
