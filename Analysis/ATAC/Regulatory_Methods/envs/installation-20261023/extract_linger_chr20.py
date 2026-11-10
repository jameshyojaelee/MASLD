"""Extract the exact15 files needed by the bounded chr20 LINGER1.110 stage."""
import hashlib,json,pathlib,tarfile,time
root=pathlib.Path('Analysis/ATAC/Regulatory_Methods/inputs/linger_1.110_zenodo23399091')
manifest=json.loads((root/'archive_manifest.json').read_text())
expected=['bulk_gene_all.txt','TFName.txt','Match2.txt','all_hg19.txt','motifWeight.txt','hg38_hg19_pair.bed','RE_gene_corr_hg19.bed','MotifTarget_Matrix_chr20.txt','MotifTarget_matrix_chr20.bed','hg38_Peaks_chr20.bed','chr20_gene.txt','chr20_index.txt','chr20_index_all.txt','fisher_chr20.pt','all_models_chr20.pt']
wanted={'data_bulk/'+name for name in expected};entries={item['path']:item for item in manifest['members']}
if not wanted.issubset(entries):raise ValueError('Missing requested files:'+str(wanted-set(entries)))
output=root/'minimal_chr20';output.mkdir(exist_ok=False);rows=[];start=time.monotonic()
with tarfile.open(root/'data_bulk.tar.gz','r|gz') as archive:
 for member in archive:
  if member.name not in wanted:continue
  path=pathlib.PurePosixPath(member.name)
  if not member.isfile() or path.is_absolute() or '..' in path.parts:raise ValueError('Unsafe requested member')
  if member.size!=entries[member.name]['bytes']:raise ValueError('Manifest/member size mismatch')
  destination=output.joinpath(*path.parts);destination.parent.mkdir(parents=True,exist_ok=True)
  if not destination.resolve().is_relative_to(output.resolve()):raise ValueError('Destination escaped subset directory')
  sha=hashlib.sha256();count=0
  with archive.extractfile(member) as source,destination.open('xb') as handle:
   for block in iter(lambda:source.read(1024*1024),b''):
    handle.write(block);sha.update(block);count+=len(block)
  if count!=member.size:raise ValueError('Truncated member')
  rows.append({'path':member.name,'bytes':count,'sha256':sha.hexdigest()});print(json.dumps(rows[-1]),flush=True)
  if len(rows)==len(wanted):break
if set(row['path'] for row in rows)!=wanted:raise ValueError('Incomplete subset extraction')
(output/'extraction.json').write_text(json.dumps({'archive_sha256':manifest['archive_sha256'],'files':rows,'extracted_bytes':sum(row['bytes'] for row in rows),'elapsed_seconds':time.monotonic()-start,'scope':'global human prior plus chromosome20 only; bounded DPM1 diagnostic','model_serialization_and_index_compatibility':'pending actual native preflight'},indent=2)+'\n')
print('EXTRACTED',len(rows),sum(row['bytes'] for row in rows),flush=True)
