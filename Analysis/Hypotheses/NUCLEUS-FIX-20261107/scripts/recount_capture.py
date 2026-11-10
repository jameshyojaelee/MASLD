"""Joint recount of a verified experiment; never combine different captures."""
import json
import os
from pathlib import Path
import subprocess
import sys

TASK = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/NUCLEUS-FIX-20261107')
assert os.environ.get('SLURM_JOB_ID')
experiment = sys.argv[1]
record, = [r for r in json.loads((TASK / 'checks/additional_recounts.json').read_text()) if r['experiment'] == experiment]
assert record['dataset'] == 'GSE202379', 'Liver Atlas paired runs include ADT and must not be jointly counted as RNA'
out = TASK / 'recovery' / experiment
out.mkdir(parents=True, exist_ok=False)
fastqs = out / 'fastqs'; fastqs.mkdir()
for row in record['fastqs']:
    (fastqs / f"{row['run']}_S1_L001_R{row['mate']}_001.fastq.gz").symlink_to(row['path'])
runs = sorted({r['run'] for r in record['fastqs']})
cmd = ['/gpfs/commons/home/jameslee/cellranger-9.0.0/cellranger', 'count',
       '--id=' + experiment, '--transcriptome=/gpfs/commons/home/jameslee/reference_genome/cellranger/refdata-gex-GRCh38-2024-A',
       '--fastqs=' + str(fastqs), '--sample=' + ','.join(runs),
       '--create-bam=false', '--nosecondary', '--chemistry=auto',
       '--localcores=8', '--localmem=56', '--output-dir=' + str(out / 'count')]
(out / 'source_and_command.json').write_text(json.dumps(dict(source=record, command=cmd), indent=2))
subprocess.run(cmd, check=True, cwd=out)
(out / 'completion.json').write_text(json.dumps(dict(library=experiment, status='completed', n_runs=len(runs)), indent=2))
