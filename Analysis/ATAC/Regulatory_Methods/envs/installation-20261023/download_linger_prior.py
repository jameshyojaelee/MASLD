"""Download one official, version-qualified LINGER1.110 resource archive."""
import datetime,hashlib,json,pathlib,time,urllib.request
root=pathlib.Path('Analysis/ATAC/Regulatory_Methods/inputs/linger_1.110_zenodo23399091')
root.mkdir(parents=True,exist_ok=False)
url='https://zenodo.org/api/records/23399091/files/data_bulk.tar.gz/content'
expected_bytes=20812511726
expected_sha='acccd10b1c8d8755da7de85adec2431052d2ceb1b8ba3adc0309b26a3fb67452'
partial=root/'data_bulk.tar.gz.partial'
sha=hashlib.sha256();size=0;start=time.monotonic();next_report=0
req=urllib.request.Request(url,headers={'Cache-Control':'no-cache'})
with urllib.request.urlopen(req,timeout=120) as response,partial.open('xb') as handle:
 if response.status!=200:raise ValueError(f'Expected complete response200, got{response.status}')
 for piece in iter(lambda:response.read(4*1024*1024),b''):
  handle.write(piece);sha.update(piece);size+=len(piece)
  if size>expected_bytes:raise ValueError('Archive larger than published metadata')
  if time.monotonic()>=next_report:
   print(json.dumps({'bytes':size,'total':expected_bytes,'elapsed_seconds':round(time.monotonic()-start,1)}),flush=True);next_report=time.monotonic()+30
if size!=expected_bytes or sha.hexdigest()!=expected_sha:raise ValueError('Published size/SHA256 mismatch')
partial.rename(root/'data_bulk.tar.gz')
(root/'download.json').write_text(json.dumps({'url':url,'record':'23399091','declared_method_version':'1.110','bytes':size,'sha256':sha.hexdigest(),'integrity':'verified','elapsed_seconds':time.monotonic()-start,'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2)+'\n')
print('VERIFIED',size,sha.hexdigest(),flush=True)
