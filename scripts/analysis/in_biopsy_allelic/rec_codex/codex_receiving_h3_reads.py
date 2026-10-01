"""Fixed seventeen-person single-read ChIP transport; no fitting or accuracy.

All RNA predictions must be complete and hash-bound before any receiving read.
Native title tokens are metadata, never inferred biological library ownership.
"""
import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time
import traceback
from urllib.parse import urlsplit
from urllib.request import urlopen

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
REC = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
ROSTER = REC / 'paired_h3k27ac_raw_roster_21997742/raw_files.tsv'
RUNS = ROSTER.with_name('raw_runs.tsv')
SAF = REC / 'fixed_h3_saf_21998535/fixed_h3_regions.saf'
INDEX = REC / 'codex_h3_genome_index_21998702'
ENV = Path('/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin')
DONORS = ('B1','B7','B8','B9','B10','B12','B15','B19','B21','B22','B24','B26','B36','B38','B41','B46','B47')
GROUPS = (('B7','B8','B9','B10'), ('B12','B15','B19','B21'),
          ('B22','B24','B26','B36'), ('B38','B41','B46','B47'))
RNA_JOBS = (21998814,21998815,21998816,21998817)
GUARDS = {ROSTER:'f91c00fb37387020a624a02a388ff00d440efc7fa5bdac5ef41da2a8b1ae7e78',
          RUNS:'69ede048c0975d464c71703193fb279bd77c9e018bf40fcbcdcd6c1ea4fe483b',
          INDEX/'summary.json':'fa2a9cee1f6566880fcf82f8bdab54b02285b697a6871146e8606a3b9a254cc6',
          INDEX/'protocol.json':'5700a859f7224d4ec6560aa576eb27a9515c458bbd05432caf14cdcd3a898efc',
          SAF:'3c0ae95792c0c67b022bddf5a702b5ac84ba04040be24030edc9b22029e3d598',
          SAF.parent/'summary.json':'feb1a674d911725ebf344438597881f06a6c4f2e08bbacbd2383cec487130237',
          HERE/'codex_build_h3_genome_index.py':'c8516a6f3beb38bfdd12f7534958fecd92fb229956e3ee3a62ede95313774a3b',
          REC/'codex_native_rna_pilot_21998663/summary.json':'97d80729969fb1f075d9a18c4134a8d9a94cd570d6f7643dbcfac6ef25c41b46'}
CAP = 250_000_000_000
STOP = 230_000_000_000
SEED = 20260930


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path, algorithm='sha256'):
    h = hashlib.new(algorithm)
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1<<20), b''):
            h.update(block)
    return h.hexdigest()


def guard(values):
    for p,h in values.items():
        require(sha(p)==h, 'Frozen input changed: '+str(p))


def read_json(path):
    require(path.is_file() and path.stat().st_size<5_000_000, 'Missing/bounded completion receipt: '+str(path))
    return json.loads(path.read_text())


def prerequisites():
    guard(GUARDS)
    frozen = dict(GUARDS)
    index = read_json(INDEX/'summary.json')
    require(index.get('index_built') and index.get('reconstructed_canonical_sequence_agrees')
            and index.get('contigs')==194 and index.get('bases')==3099750718, 'Genome index not fully verified')
    protocol = read_json(INDEX/'protocol.json')
    require(sha(INDEX/'protocol.json')==index['protocol_sha256']
            and protocol['source_sha256']==GUARDS[HERE/'codex_build_h3_genome_index.py'], 'Genome index lineage differs')
    require(index['index_prefix']==str(INDEX/'GRCh38_ensembl98_primary') and len(index['index_files'])==6,
            'Wrong index path or incomplete index')
    for p,info in index['index_files'].items():
        p=Path(p)
        require(p.parent==INDEX and p.stat().st_size==info['bytes'], 'Index escaped/missing payload')
        frozen[p]=info['sha256']
    frozen[INDEX/'summary.json']=sha(INDEX/'summary.json')
    frozen[INDEX/'protocol.json']=index['protocol_sha256']
    summaries={'B1':REC/'codex_native_rna_pilot_21998663/summary.json'}
    for job,group in zip(RNA_JOBS,GROUPS):
        base=REC/('codex_native_rna_batch_'+str(job))
        result=read_json(base/'summary.json')
        require(result.get('group_complete') and tuple(result['donors'])==group
                and set(result['status'])==set(group), 'RNA group not complete/exact')
        frozen[base/'summary.json']=sha(base/'summary.json')
        for donor in group:
            status=result['status'][donor]
            require(status['state']=='succeeded' and Path(status['output'])==base/donor, 'RNA donor not complete/exact')
            summaries[donor]=base/donor/'summary.json'
            frozen[summaries[donor]]=status['summary_sha256']
    require(set(summaries)==set(DONORS), 'Require all seventeen fixed RNA donors')
    for donor,p in summaries.items():
        s=read_json(p)
        require(s['donor']==donor and s.get('native_quantification_completed') and s.get('frozen_predictions_completed')
                and s.get('modeled_genes')==42163 and s.get('native_genes')==60623 and not s.get('rows_fabricated')
                and not s.get('receiving_H3_read') and not s.get('model_fitted') and not s.get('predictive_accuracy_evaluated'),
                'RNA completion/measurement boundary differs')
        frozen[p]=frozen.get(p,sha(p))
        frozen[p.parent/'protocol.json']=s['protocol_sha256']
        require(set(s['predictions'])=={'global','local'}, 'Missing paired model predictions')
        for name,form in [('global','rrr'),('local','rrr_offset_cis')]:
            prediction=s['predictions'][name]; payload=Path(prediction['path'])
            require(payload.parent==p.parent and prediction['form']==form and prediction['transport']=='raw',
                    'RNA prediction form/transport/path differs')
            frozen[payload]=prediction['sha256']
    # Hash payload bytes only, never deserialize RNA or prediction arrays.
    guard(frozen)
    with ROSTER.open(newline='') as f:
        all_rows=list(csv.DictReader(f,delimiter='\t'))
    rows=[r for r in all_rows if r['assay']=='H3K27ac']
    require(set(r['donor'] for r in rows)==set(DONORS) and len(rows)==202
            and len({r['run_accession'] for r in rows})==202
            and len({r['fastq_url'] for r in rows})==202
            and all(r['file_number']=='1' for r in rows)
            and sum(int(r['fastq_bytes']) for r in rows)==56783483922,
            'Exact H3-only single-file roster differs')
    with RUNS.open(newline='') as f:
        native_runs=[r for r in csv.DictReader(f,delimiter='\t') if r['assay']=='H3K27ac']
    require(len(native_runs)==202 and len({r['run_accession'] for r in native_runs})==202,
            'Missing or duplicate native H3 run metadata')
    by_run={r['run_accession']:r for r in native_runs}
    require(set(by_run)=={r['run_accession'] for r in rows}, 'H3 file/run identity axes differ')
    public_fields=('run_accession','donor','GSM','experiment_accession','source_biosample')
    for row in rows:
        native=by_run[row['run_accession']]
        require(all(native[k]==row[k] and row[k] for k in public_fields)
                and native['sample_accession']==row['source_biosample'], 'H3 run/source ownership join differs')
        require(native['library_strategy']=='ChIP-Seq' and native['library_layout']=='SINGLE'
                and native['fastq_file_count']=='1', 'Require explicit SINGLE ChIP-Seq one-file layout')
        require(native['fastq_ftp']==row['fastq_url'] and native['fastq_md5']==row['fastq_md5']
                and int(native['fastq_bytes'])==int(row['fastq_bytes'])==int(native['total_fastq_bytes']),
                'Native H3 one-file URL/size/checksum join differs')
    return frozen,rows,index


def output_bytes(path):
    return sum(p.stat().st_size for p in path.rglob('*') if p.is_file())


def bound(output, donor_output, donor_cap):
    require(output_bytes(output)<=STOP and output_bytes(donor_output)<=donor_cap,
            'Sampled saved-output stop exceeded; partials preserved')


def write_json(path,value):
    with path.open('x') as f:
        json.dump(value,f,indent=2);f.write('\n')


def acquire(row,path,output,donor_output,donor_cap):
    url='https://'+row['fastq_url'];parsed=urlsplit(url)
    require(parsed.hostname=='ftp.sra.ebi.ac.uk' and parsed.path.startswith('/vol1/fastq/')
            and not parsed.query and not parsed.fragment, 'Unexpected read URL')
    expected=int(row['fastq_bytes']); size=0; md5=hashlib.md5(); digest=hashlib.sha256()
    partial=Path(str(path)+'.partial')
    with urlopen(url,timeout=120) as response,partial.open('xb') as dest:
        require(response.geturl()==url, 'Unexpected read redirect')
        length=response.headers.get('Content-Length')
        require(length is not None and int(length)==expected, 'Missing/different frozen Content-Length')
        for block in iter(lambda:response.read(1<<20),b''):
            size+=len(block);require(size<=expected, 'Read file exceeds frozen size')
            dest.write(block);md5.update(block);digest.update(block)
            if size % (64 << 20) < len(block):
                require(output_bytes(output)<=CAP and output_bytes(donor_output)<=donor_cap, 'Download cap exceeded; partial retained')
    require(output_bytes(output)<=CAP and output_bytes(donor_output)<=donor_cap, 'Download cap exceeded; partial retained')
    require(size==expected and md5.hexdigest()==row['fastq_md5'], 'Read size/MD5 differs; partial retained')
    expanded=0
    with gzip.open(partial,'rb') as f:
        for block in iter(lambda:f.read(1<<20),b''):
            expanded+=len(block);require(expanded<=50*expected, 'Read decompression cap exceeded')
    require(not path.exists(), 'Refuse existing read overwrite');partial.rename(path)
    return dict(**row,path=str(path),sha256=digest.hexdigest(),expanded_bytes=expanded,gzip_crc_verified=True)


def run(commands,log,output,donor_output,donor_cap):
    """Pipe alignment directly to BAM; check every process and preserve failures."""
    children=[]
    with log.open('x') as handle:
        try:
            previous=None
            for i,cmd in enumerate(commands):
                child=subprocess.Popen(cmd,stdin=previous,stdout=subprocess.PIPE if i<len(commands)-1 else handle,
                                       stderr=handle,start_new_session=True)
                if previous is not None:previous.close()
                previous=child.stdout;children.append(child)
            while any(c.poll() is None for c in children):
                bound(output,donor_output,donor_cap)
                require(not any(c.poll() not in (None,0) for c in children), 'Pipeline failed: '+str(log))
                time.sleep(5)
            require(all(c.wait()==0 for c in children), 'Command failed: '+str(log))
            bound(output,donor_output,donor_cap)
        finally:
            for child in children:
                if child.poll() is None:os.killpg(child.pid,signal.SIGTERM)
            for child in children:
                try:child.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    os.killpg(child.pid,signal.SIGKILL);child.wait()


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.parse_args()
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and int(os.environ.get('SLURM_CPUS_PER_TASK','0'))>=8,
            'Receiving read processing requires eight allocated CPUs')
    frozen,rows,index=prerequisites() # Fail before any H3 download/read if RNA/index incomplete.
    bowtie=shutil.which('bowtie2');samtools=ENV/'samtools';counts=ENV/'featureCounts'
    require(bowtie, 'Load Bowtie2/2.5.4-linux-x86_64')
    versions={}
    for name,cmd in [('bowtie2',[bowtie,'--version']),('samtools',[str(samtools),'--version']),('featureCounts',[str(counts),'-v'])]:
        p=subprocess.run(cmd,capture_output=True,text=True,check=True);versions[name]=p.stdout+p.stderr
    require('version 2.5.4' in versions['bowtie2'] and 'samtools 1.22.1' in versions['samtools']
            and '2.1.1' in versions['featureCounts'], 'Installed tool versions changed')
    frozen.update({Path(__file__):sha(Path(__file__)), HERE/'run_codex_receiving_h3_reads.sbatch':sha(HERE/'run_codex_receiving_h3_reads.sbatch'),
                   Path(bowtie):sha(bowtie),samtools:sha(samtools),counts:sha(counts)})
    output=REC/('codex_receiving_h3_reads_'+os.environ['SLURM_JOB_ID']);output.mkdir(exist_ok=False)
    donor_caps={d:5*sum(int(r['fastq_bytes']) for r in rows if r['donor']==d)+2_000_000_000 for d in DONORS}
    write_json(output/'protocol.json',dict(fixed_donors=DONORS,h3_only_files=202,compressed_bytes=56783483922,
        input_sha256={str(p):h for p,h in frozen.items()},versions=versions,donor_caps=donor_caps,
        sampled_stop_bytes=STOP,hard_cap_bytes=CAP,cap_monitor_seconds=5,alignment_seed=SEED,
        alignment='Bowtie2 sensitive end-to-end; Phred33; single unpaired default best primary alignment; no trimming or library pooling/dedup',
        count_unit='mapped primary nonsupplementary single reads MAPQ>=10; unstranded one-region assignment; no -O/-M/Input/PCRdedup/peakcalls',
        interpretation='Qualified native single-read ChIP transport from CUT&RUN-trained residual predictor; not original fragment/PCR/private-library equivalence',
        RNA_predictions_frozen_before_H3=True,model_fit=False,predictive_accuracy_evaluated=False))
    (output/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    completed=[];donor=None;phase='initialization'
    try:
        for donor in DONORS:
            guard(frozen);directory=output/donor;directory.mkdir();raw=directory/'raw';raw.mkdir()
            bindings={};filters={};acquired=[];commands=[]
            for row in sorted((r for r in rows if r['donor']==donor),key=lambda r:r['run_accession']):
                phase='acquisition';run_id=row['run_accession'];fastq=raw/Path(row['fastq_url']).name
                acquired.append(acquire(row,fastq,output,directory,donor_caps[donor]))
                write_json(directory/(run_id+'_acquisition.json'),acquired[-1])
                aligned=directory/(run_id+'.aligned.bam');filtered=directory/(run_id+'.filtered.bam')
                align=[bowtie,'--end-to-end','--sensitive','--phred33','--seed',str(SEED),'-p','7',
                       '-x',index['index_prefix'],'-U',str(fastq),'--rg-id',run_id,'--rg','SM:'+donor]
                to_bam=[str(samtools),'view','-b','-o',str(aligned),'-']
                phase='alignment';run([align,to_bam],directory/(run_id+'_alignment.log'),output,directory,donor_caps[donor])
                prefilter=[str(samtools),'view','-b','-F','2308','-q','10','-o',str(filtered),str(aligned)]
                phase='primary_filter';run([prefilter],directory/(run_id+'_filter.log'),output,directory,donor_caps[donor])
                run([[str(samtools),'quickcheck',str(aligned),str(filtered)]],directory/(run_id+'_quickcheck.log'),output,directory,donor_caps[donor])
                bindings[str(filtered)]={k:row[k] for k in ('run_accession','donor','GSM','experiment_accession','source_biosample')}
                filters[str(filtered)]=prefilter;commands.append(dict(alignment=align,to_bam=to_bam,prefilter=prefilter))
            phase='fixed_region_counting';count_path=directory/'native_featurecounts.tsv'
            command=[str(counts),'-T','8','-F','SAF','-a',str(SAF),'-s','0','-Q','10','--primary','-o',str(count_path),*bindings]
            run([command],directory/'featurecounts.log',output,directory,donor_caps[donor]);guard(frozen)
            write_json(directory/'summary.json',dict(donor=donor,h3_reads_processed=True,count_table=str(count_path),
                count_command=command,columns=bindings,prefilter_commands=filters,alignment_commands=commands,
                count_sha256=sha(count_path),count_summary_sha256=sha(Path(str(count_path)+'.summary')),
                native_read_files=acquired,saved_output_bytes=output_bytes(directory),predictions_refitted=False,
                count_axis_assembly_or_target_evaluation_run=False))
            completed.append(donor);print(json.dumps(dict(donor=donor,completed=len(completed),total=17)),flush=True)
        guard(frozen)
        write_json(output/'summary.json',dict(completed=DONORS,all17_reads_processed=True,biological_n=17,
            saved_output_bytes=output_bytes(output),protocol_sha256=sha(output/'protocol.json'),
            donor_summary_sha256={d:sha(output/d/'summary.json') for d in DONORS},predictive_accuracy_evaluated=False))
    except BaseException:
        write_json(output/'failure.json',dict(donor=donor,phase=phase,completed=completed,partials_preserved=True,traceback=traceback.format_exc()))
        raise


if __name__=='__main__':
    main()
