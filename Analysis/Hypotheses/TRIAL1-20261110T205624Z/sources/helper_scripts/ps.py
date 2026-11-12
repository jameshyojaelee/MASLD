import json,sys,urllib.request,urllib.parse,time
q=sys.argv[1]; n=sys.argv[2] if len(sys.argv)>2 else '10'
u='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmode=json&retmax=%s&term=%s'%(n,urllib.parse.quote(q))
ids=json.load(urllib.request.urlopen(u))['esearchresult']['idlist']
if not ids: print('none'); sys.exit()
time.sleep(0.4)
d=json.load(urllib.request.urlopen('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id='+','.join(ids)))['result']
for k in d['uids']: print(k, d[k]['pubdate'][:8], d[k]['source'][:20],'|', d[k]['title'][:170])
