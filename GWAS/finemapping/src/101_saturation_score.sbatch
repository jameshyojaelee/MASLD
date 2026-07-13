#!/usr/bin/env bash
#SBATCH --job-name=chromBPNet
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --array=0-26%4
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/101_saturation_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/101_saturation_%A_%a.err
# Score saturation SNVs with an independently trained, open-weight ChromBPNet model.
# Required submission exports: MODEL, PEAKS, CHROMSIZES. This script intentionally
# has no AlphaGenome path and cannot modify canonical PIPs/convergence.
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SEQ=$ROOT/GWAS/finemapping/results/seqfunc/haplotype_saturation/sequences
VS=$ROOT/GWAS/finemapping/src/chrombpnet_variant_scorer
: "${MODEL:?export MODEL=/absolute/path/model.h5}"
: "${PEAKS:?export PEAKS=/absolute/path/peaks.bed}"
: "${CHROMSIZES:?export CHROMSIZES=/absolute/path/chrom.sizes}"
: "${FOLD_JSON:?export FOLD_JSON=/absolute/path/fold.json}"
: "${MODEL_FOLD_ID:?export MODEL_FOLD_ID=foldN}"
SHARD=$(printf '%03d' "$SLURM_ARRAY_TASK_ID")
LIST=$SEQ/genomic_variants.shard${SHARD}.chrombpnet.tsv
OUT=$SEQ/scores; mkdir -p "$OUT"
for f in "$MODEL" "$PEAKS" "$CHROMSIZES" "$FOLD_JSON" "$LIST"; do test -s "$f" || { echo "FATAL missing $f"; exit 2; }; done
ENV=/gpfs/commons/home/jameslee/micromamba/envs/chrombpnet_vs
export PYTHONNOUSERSITE=1 TF_FORCE_GPU_ALLOW_GROWTH=true
# A chromosome may only be scored by a model for which it was a held-out test
# chromosome. This prevents optimistic in-sample saturation maps.
"$ENV/bin/python" - "$LIST" "$FOLD_JSON" <<'PY'
import csv,json,sys
query={r[0] for r in csv.reader(open(sys.argv[1]),delimiter="\t") if r}
fold=json.load(open(sys.argv[2])); held=set(fold.get("test",[]))
bad=sorted(query-held)
if bad: raise SystemExit("FATAL training/validation chromosome scoring prohibited; not in fold test set: "+",".join(bad))
print("held-out chromosome gate PASS:",sorted(query))
PY
"$ENV/bin/python" "$VS/src/variant_scoring.py" --list "$LIST" \
  --genome /gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A/fasta/genome.fa \
  --model "$MODEL" --out_prefix "$OUT/genomic_shard${SHARD}" --chrom_sizes "$CHROMSIZES" \
  --peaks "$PEAKS" --max_peaks 30000 --total_shuf 30000 --batch_size 256 \
  --schema chrombpnet --random_seed 1234 --no_hdf5
test -s "$OUT/genomic_shard${SHARD}.variant_scores.tsv"
"$ENV/bin/python" - "$MODEL" "$FOLD_JSON" "$OUT/genomic_shard${SHARD}.provenance.json" "$MODEL_FOLD_ID" <<'PY'
import hashlib,json,sys
def sha(p):
 h=hashlib.sha256()
 with open(p,"rb") as f:
  for b in iter(lambda:f.read(8<<20),b""): h.update(b)
 return h.hexdigest()
json.dump({"model":sys.argv[1],"model_sha256":sha(sys.argv[1]),"fold_json":sys.argv[2],"model_fold_id":sys.argv[4],
 "fold_json_sha256":sha(sys.argv[2]),"heldout_test_chromosomes":json.load(open(sys.argv[2]))["test"],
 "interpretation":"mechanistic_nomination_not_causal_evidence","apply_only_firewall":True},
 open(sys.argv[3],"w"),indent=2)
PY
