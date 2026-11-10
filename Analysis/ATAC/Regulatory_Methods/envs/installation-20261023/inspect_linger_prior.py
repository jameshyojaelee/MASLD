"""Read complete verified archive member metadata without extracting any files."""
import hashlib,json,pathlib,tarfile,time
root=pathlib.Path('Analysis/ATAC/Regulatory_Methods/inputs/linger_1.110_zenodo23399091')
receipt=json.loads((root/'download.json').read_text())
assert receipt['integrity']=='verified'
manifest=root/'archive_manifest.json'
if manifest.exists():raise FileExistsError(manifest)
start=time.monotonic();members=[]
with tarfile.open(root/'data_bulk.tar.gz','r|gz') as archive:
 for member in archive:
  path=pathlib.PurePosixPath(member.name)
  if path.is_absolute() or '..' in path.parts or not path.parts or path.parts[0]!='data_bulk':raise ValueError('Unsafe member path:'+member.name)
  if member.issym() or member.islnk() or not (member.isfile() or member.isdir()):raise ValueError('Unexpected member type:'+member.name)
  members.append({'path':member.name,'bytes':member.size,'regular_file':member.isfile()})
  print(json.dumps(members[-1]),flush=True)
result={'archive_sha256':receipt['sha256'],'path_validation':'all members within data_bulk, no links or special files','members':members,'n_members':len(members),'expanded_regular_file_bytes':sum(x['bytes'] for x in members if x['regular_file']),'elapsed_seconds':time.monotonic()-start,'extraction_run':False}
manifest.write_text(json.dumps(result,indent=2)+'\n')
print('MANIFEST_SUMMARY',result['n_members'],result['expanded_regular_file_bytes'],result['elapsed_seconds'],flush=True)
