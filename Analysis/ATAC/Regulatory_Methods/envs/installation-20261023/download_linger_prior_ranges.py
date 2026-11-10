"""Four-stream resumable official prior download; no extraction or model loading."""
import concurrent.futures,datetime,hashlib,json,os,pathlib,threading,time,urllib.request,urllib.error,email.utils
root=pathlib.Path('Analysis/ATAC/Regulatory_Methods/inputs/linger_1.110_zenodo23399091')
partial=root/'data_bulk.tar.gz.partial';final=root/'data_bulk.tar.gz';progress=root/'range_progress.json'
url='https://zenodo.org/api/records/23399091/files/data_bulk.tar.gz/content'
total=20812511726;expected='acccd10b1c8d8755da7de85adec2431052d2ceb1b8ba3adc0309b26a3fb67452'
if final.exists():raise FileExistsError(final)
state=json.loads(progress.read_text()) if progress.exists() else {'verified_initial_sequential_prefix_bytes':partial.stat().st_size,'completed_ranges':[]}
start=state['verified_initial_sequential_prefix_bytes'];block=32*1024*1024
ranges=[(a,min(a+block-1,total-1)) for a in range(start,total,block)]
complete={tuple(x) for x in state['completed_ranges']};lock=threading.Lock();clock=time.monotonic();fd=os.open(partial,os.O_RDWR)
def save():
 temporary=progress.with_suffix('.tmp');temporary.write_text(json.dumps(state)+'\n');temporary.replace(progress)
save()
def download(interval):
 if interval in complete:return
 a,b=interval
 for attempt in range(8):
  try:
   request=urllib.request.Request(url,headers={'Range':f'bytes={a}-{b}','Cache-Control':'no-cache'})
   with urllib.request.urlopen(request,timeout=120) as response:
    if response.status!=206 or response.headers.get('Content-Range')!=f'bytes {a}-{b}/{total}':raise ValueError('Content-Range mismatch')
    offset=a
    while offset<=b:
     data=response.read(min(1024*1024,b-offset+1))
     if not data:raise OSError('truncated range')
     written=os.pwrite(fd,data,offset)
     if written!=len(data):raise OSError('short disk write')
     offset+=written
   with lock:
    state['completed_ranges'].append([a,b]);save()
    done=start+sum(y-x+1 for x,y in state['completed_ranges'])
    print(json.dumps({'bytes':done,'total':total,'elapsed_seconds':round(time.monotonic()-clock,1),'streams':4}),flush=True)
   return
  except urllib.error.HTTPError as exc:
   if exc.code not in [429,500,502,503,504]:raise
   retry=exc.headers.get('Retry-After','60')
   try:delay=float(retry)
   except ValueError:delay=max(0,email.utils.parsedate_to_datetime(retry).timestamp()-time.time())
   print(f'HTTP{exc.code}; honoring Retry-After={delay}',flush=True);time.sleep(max(delay,2**attempt))
  except (TimeoutError,OSError):
   if attempt==7:raise
   time.sleep(min(60,2**attempt))
 raise RuntimeError('range retries exhausted')
try:
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:list(pool.map(download,ranges))
 os.fsync(fd)
finally:os.close(fd)
if partial.stat().st_size!=total:raise ValueError('size mismatch')
sha=hashlib.sha256()
with partial.open('rb') as handle:
 for data in iter(lambda:handle.read(8*1024*1024),b''):sha.update(data)
if sha.hexdigest()!=expected:raise ValueError('published SHA256 mismatch')
partial.rename(final)
(root/'download.json').write_text(json.dumps({'url':url,'record':'23399091','declared_method_version':'1.110','bytes':total,'sha256':sha.hexdigest(),'integrity':'verified','streams':4,'elapsed_seconds_range_phase':time.monotonic()-clock,'resumed_prefix_bytes':start,'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2)+'\n')
print('VERIFIED',sha.hexdigest(),flush=True)
