#!/usr/bin/env bash
#SBATCH --job-name=chromBPNet
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=90:00:00
# Relative, resolved against the submission working directory, which
# submit_saturation_v2.sh sets to $SATURATION_V2_ROOT.  These were absolute
# paths into haplotype_saturation/v2/logs/, so a v3 saturation run would have
# written its logs into the v2 tree that the frozen release still references.
#SBATCH --output=logs/score_%A_%a.out
#SBATCH --error=logs/score_%A_%a.err
set -euo pipefail
ROOT=${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}
OUT=${SATURATION_V2_ROOT:-$ROOT/GWAS/finemapping/results/seqfunc/haplotype_saturation/v2}
SEQ=$OUT/sequences; SCORES=$OUT/scores; mkdir -p "$SCORES" "$OUT/logs"
GATE=$OUT/gates/upstream_gate.json; MODELS=$OUT/gates/model_manifest.locked.tsv; SHARDS=$SEQ/scoring_shards.tsv
for f in "$GATE" "$MODELS" "$SHARDS"; do test -s "$f" || { echo "FATAL missing $f"; exit 2; }; done
python3 - "$GATE" <<'PY'
import json,sys
x=json.load(open(sys.argv[1]))
if not x.get("saturation_submission_allowed",False): raise SystemExit("upstream gate BLOCKED")
PY
ROW=$(awk -F '\t' -v i="$SLURM_ARRAY_TASK_ID" 'NR==i+2{print;exit}' "$SHARDS" | tr -d '\r')
test -n "$ROW" || { echo "FATAL no shard row $SLURM_ARRAY_TASK_ID"; exit 2; }
IFS=$'\t' read -r IDX SHARD FOLD TEST_CHROMS NVAR TABLE LIST <<< "$ROW"
MROW=$(awk -F '\t' -v f="$FOLD" 'NR>1&&$1==f{print;exit}' "$MODELS" | tr -d '\r')
test -n "$MROW" || { echo "FATAL no locked model for $FOLD"; exit 2; }
IFS=$'\t' read -r MFOLD MTEST MODEL MODEL_SHA PEAKS PEAKS_SHA FOLD_JSON FOLD_SHA REFERENCE_FASTA LOCKED_CHROMSIZES LOCKED_SIZES_SHA <<< "$MROW"
for f in "$MODEL" "$PEAKS" "$FOLD_JSON" "$LIST"; do test -s "$f" || { echo "FATAL missing $f"; exit 2; }; done
python3 - "$LIST" "$FOLD_JSON" "$FOLD" <<'PY'
import csv,json,sys
chroms={r[0] for r in csv.reader(open(sys.argv[1]),delimiter="\t") if r}
f=json.load(open(sys.argv[2])); held=set(f.get("test",f.get("test_chromosomes",[])))
if not chroms or not chroms<=held: raise SystemExit(f"heldout routing failure {chroms-held}")
print(sys.argv[3],sorted(chroms),"heldout PASS")
PY
ENV=${CHROMBPNET_SCORER_ENV:-/gpfs/commons/home/jameslee/micromamba/envs/chrombpnet_vs}
VS=${CHROMBPNET_VARIANT_SCORER:-$ROOT/GWAS/finemapping/src/chrombpnet_variant_scorer}
FASTA=${GRCH38_FASTA:-$REFERENCE_FASTA}
CHROMSIZES=${GRCH38_CHROMSIZES:-$LOCKED_CHROMSIZES}
export PYTHONNOUSERSITE=1 TF_FORCE_GPU_ALLOW_GROWTH=true
PREFIX=$SCORES/$SHARD
"$ENV/bin/python" "$VS/src/variant_scoring.py" --list "$LIST" --genome "$FASTA" \
  --model "$MODEL" --out_prefix "$PREFIX" --chrom_sizes "$CHROMSIZES" --peaks "$PEAKS" \
  --max_peaks 30000 --total_shuf 30000 --batch_size 256 --schema chrombpnet --random_seed 1234 --no_hdf5
test -s "$PREFIX.variant_scores.tsv"
python3 - "$PREFIX.provenance.json" "$SHARD" "$FOLD" "$MODEL" "$MODEL_SHA" "$FOLD_JSON" "$FOLD_SHA" "$TEST_CHROMS" "$NVAR" "$FASTA" "$CHROMSIZES" "$LOCKED_SIZES_SHA" <<'PY'
import json,os,sys
out,shard,fold,model,msha,fjson,fsha,chroms,n,fasta,sizes,sizes_sha=sys.argv[1:]
tmp=out+f".tmp.{os.getpid()}"
json.dump({"status":"PASS","shard_id":shard,"fold_id":fold,"heldout_test_chromosomes":chroms.split(";"),
 "n_expected":int(n),"model":model,"model_sha256":msha,"fold_json":fjson,"fold_json_sha256":fsha,
 "reference_fasta":fasta,"chrom_sizes":sizes,"chrom_sizes_sha256":sizes_sha,
 "interpretation":"mechanistic_nomination_not_causal_evidence","apply_only_firewall":True},open(tmp,"w"),indent=2)
os.replace(tmp,out)
PY
