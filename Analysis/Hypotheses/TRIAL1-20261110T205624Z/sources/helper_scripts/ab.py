import sys,urllib.request,xml.etree.ElementTree as ET
ids=sys.argv[1]
x=urllib.request.urlopen('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&retmode=xml&id='+ids).read()
open(sys.argv[2] if len(sys.argv)>2 else '/dev/null','ab').write(x)
r=ET.fromstring(x)
for a in r.findall('.//PubmedArticle'):
    pmid=a.findtext('.//PMID'); t=''.join(a.find('.//ArticleTitle').itertext())
    doi=[e.text for e in a.findall('.//ArticleId') if e.get('IdType')=='doi']
    j=a.findtext('.//ISOAbbreviation'); y=a.findtext('.//PubDate/Year') or a.findtext('.//PubDate/MedlineDate')
    print('##',pmid,doi[:1],j,y,'|',t)
    for s in a.findall('.//Abstract/AbstractText'):
        txt=''.join(s.itertext())
        if s.get('Label') in ('BACKGROUND','INTRODUCTION','BACKGROUND & AIMS','FUNDING','BACKGROUND AND AIMS'): continue
        print(' ',(s.get('Label') or '')+':',txt)
    print()
