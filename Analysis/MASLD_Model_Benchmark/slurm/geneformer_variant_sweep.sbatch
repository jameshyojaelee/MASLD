#!/usr/bin/env bash
#SBATCH --job-name=masld-gf-sweep
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --gres=gpu:l40s:1
#SBATCH --time=08:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
cd "$ROOT"
eval "$(micromamba shell hook -s bash)"
micromamba activate decima
OUT="executions/geneformer-sweep/${SLURM_JOB_ID}"
mkdir -p "$OUT"
python -c "import torch;print('cuda',torch.cuda.is_available(),torch.cuda.get_device_name(0) if torch.cuda.is_available() else '')"

BUNDLES="$ROOT/executions/geneformer-acquisition-21014827/bundles"
H5="$ROOT/executions/resource-atlas-smoke-20260821T221938Z-20984634/subset/resource_atlas_geneformer_smoke_1000.h5ad"

for MODEL in geneformer_v1_10m geneformer_v2_104m geneformer_v2_316m; do
for HEAD in linear two_layer_mlp; do
  TAG="${MODEL}__${HEAD}"
  W="$OUT/$TAG"; mkdir -p "$W"
  python - "$MODEL" "$HEAD" "$BUNDLES" "$H5" "$W" <<'PY'
import json,hashlib,pathlib,sys
model,head,bundles,h5,w = sys.argv[1:6]
b=f"{bundles}/{model}--11f311a5078fc4c8"
def art(p,role):
    d=pathlib.Path(p).read_bytes()
    return {"path":str(pathlib.Path(p).resolve()),"sha256":hashlib.sha256(d).hexdigest(),
            "size_bytes":len(d),"media_type":"application/octet-stream","role":role}
hp={"embedding_policy":"last_hidden_state_mean_over_non_special_nonpadding_tokens",
    "common_head_id":head,"outer_folds":5,
    "join_namespace":"resource_atlas_current:cell_state_mapping:donor_outer:v1",
    "cell_state_roster":["cholangiocyte","endothelial","hepatocyte","immune","mesenchymal_stromal"],
    "embedding_batch_size":32,"head_epochs":60,"device":"auto"}
base={"schema_version":"masld-bench-adapter-request-v1",
 "run_id":hashlib.sha256(f"gf-{model}-{head}".encode()).hexdigest(),
 "fit_dataset_ids":["resource_atlas_current"],"action_dataset_ids":["resource_atlas_current"],
 "run_spec":{"task_id":"cell_state_mapping","split_id":"donor_outer","stage":"smoke",
  "adaptation_regime":"common_lane","runtime_id":"gpu_geneformer","model_id":model,
  "seed":1103,"fold":0,"dataset_ids":["resource_atlas_current"],"hyperparameters":hp,
  "inputs":[art(h5,"dataset_view_data:resource_atlas_geneformer_smoke_1000_v1"),
            art(f"{b}/bundle.json",f"model_checkpoint_bundle:{model}")]}}
for a in ("prepare","fit","predict"):
    r=dict(base); r["action"]=a
    pri={"prepare":f"{w}/prepare"}
    if a=="predict": pri["fit"]=f"{w}/fit"
    if a!="prepare": r["prior_action_outputs"]=pri
    json.dump(r,open(f"{w}/req_{a}.json","w"))
PY
  echo "=== $TAG ==="
  for A in prepare fit predict; do
    PYTHONPATH="$ROOT/src/masld_bench/adapters" python \
      "$ROOT/src/masld_bench/adapters/geneformer_common_lane.py" \
      --action "$A" --request "$W/req_$A.json" --output "$W/$A" 2>&1 | grep -viE "warn|initialized|down-stream" || true
  done
  echo "DONE $TAG"
done
done
echo "SWEEP COMPLETE"
