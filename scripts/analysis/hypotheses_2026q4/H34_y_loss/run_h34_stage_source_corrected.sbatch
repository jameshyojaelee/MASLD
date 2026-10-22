#!/bin/bash
#SBATCH --job-name=H34_stage_source
#SBATCH --account=nslab
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=90:00:00
#SBATCH --output=logs/H34_stage_source_%j.out
#SBATCH --error=logs/H34_stage_source_%j.err
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
H34_SOURCE_ROOT="$PWD/Analysis/Hypotheses/H34-stage-source-corrected-20261014_${SLURM_JOB_ID}"
H34_OLD_ROOT="$PWD/Analysis/Hypotheses/H34-20261008T015507Z"
H34_CACHE_ROOT="$PWD/Analysis/Hypotheses/H34-20261008T014810Z/intermediate"
H34_ROSTER="$PWD/Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014/metadata/current_source_diagnosis_donor_roster.tsv"
H34_SCRIPT="$PWD/scripts/analysis/hypotheses_2026q4/H34_y_loss/h34_y_loss.py"
H34_PYTHON=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
mkdir "$H34_SOURCE_ROOT"
export H34_SOURCE_ROOT H34_OLD_ROOT H34_CACHE_ROOT H34_ROSTER H34_SCRIPT
export PYTHONNOUSERSITE=1 PYTHONHASHSEED=34 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cp "$H34_OLD_ROOT/PRESPEC_H34.md" "$H34_OLD_ROOT/DEVIATIONS_H34.md" "$H34_SOURCE_ROOT/"
cp -a "$H34_OLD_ROOT/geo_sample_characteristics" "$H34_SOURCE_ROOT/"
cp scripts/analysis/hypotheses_2026q4/H34_y_loss/run_h34_stage_source_corrected.sbatch "$H34_SOURCE_ROOT/"
"$H34_PYTHON" - <<'PY'
from pathlib import Path
import hashlib,json,os
out=Path(os.environ['H34_SOURCE_ROOT']);src=Path(os.environ['H34_SCRIPT']);s=src.read_text()
assert hashlib.sha256(src.read_bytes()).hexdigest()=='f7a0bfb7df1aa8a961c8073e105bcf7a2729988c54b9af2e84a60bb20f7a2f38'
a=s.index('ROSTER = REPO / (');b=s.index('\nGTF =',a)
s=s[:a]+'ROSTER = Path('+repr(os.environ['H34_ROSTER'])+')'+s[b:]
(out/'h34_y_loss_qualified_stage.py').write_text(s)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for block in iter(lambda:f.read(1048576),b''):h.update(block)
 return h.hexdigest()
paths=[src,Path(os.environ['H34_ROSTER'])]+[Path(os.environ['H34_CACHE_ROOT'])/n for n in ['pass1.npz','pass2.npz']]+[Path(os.environ['H34_OLD_ROOT'])/n for n in ['control_gene_sets.tsv','donor_celltype_excess_Yneg.tsv','donor_sex_calls.tsv','PRESPEC_H34.md']]
assert all(p.exists() for p in paths)
(out/'input_sha256_before.json').write_text(json.dumps({str(p):sha(p) for p in paths},indent=2))
PY
"$H34_PYTHON" -m pip freeze > "$H34_SOURCE_ROOT/pip_freeze.txt"
"$H34_PYTHON" "$H34_SOURCE_ROOT/h34_y_loss_qualified_stage.py" "$PWD" "$H34_SOURCE_ROOT" "$H34_CACHE_ROOT"
"$H34_PYTHON" - <<'PY'
from pathlib import Path
import hashlib,json,os,pandas as pd
from pandas.testing import assert_frame_equal
out=Path(os.environ['H34_SOURCE_ROOT']);old=Path(os.environ['H34_OLD_ROOT'])
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for block in iter(lambda:f.read(1048576),b''):h.update(block)
 return h.hexdigest()
before=json.loads((out/'input_sha256_before.json').read_text());after={p:sha(Path(p)) for p in before}
assert before==after
(out/'input_sha256_after.json').write_text(json.dumps(after,indent=2))
for n,drop in [('control_gene_sets.tsv',[]),('donor_celltype_excess_Yneg.tsv',['stage','disease_stage_coarse']),('donor_sex_calls.tsv',['stage','disease_stage_coarse']),('within_donor_summary.tsv',[])]:
 a=pd.read_csv(old/n,sep='\t').drop(columns=drop);b=pd.read_csv(out/n,sep='\t').drop(columns=drop)
 assert_frame_equal(a,b,check_exact=False,rtol=1e-10,atol=1e-10)
a=pd.read_csv(old/'stage_association.tsv',sep='\t');b=pd.read_csv(out/'stage_association.tsv',sep='\t')
assert_frame_equal(a[['cell_type','family','n_donors','tested']],b[['cell_type','family','n_donors','tested']])
a.merge(b,on='cell_type',suffixes=('_old','_corrected')).to_csv(out/'stage_estimate_comparison.tsv',sep='\t',index=False)
(out/'CHECKED.txt').write_text('Native H34 fixed model/family refit with qualified source-stage labels; unchanged source/cache hashes, donor Y-loss scores, sex calls,50 control gene sets, within-donor summary and donor counts checked. Candidate only.\n')
print(b[['cell_type','n_donors','beta_pp_per_stage','p','q_BH_immune','positive_call']].to_string(index=False))
PY
