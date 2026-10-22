#!/bin/bash
#SBATCH --job-name=H27_kidney_hf_coloc
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=180G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H27-20261008T230909Z/logs/slurm-%j.out

# H27: eGFR, CKD and heart failure colocalization at anchored MASLD loci.
# Implements PRESPEC_H27.md. One job: download outcomes, prep, 4 parallel coloc
# workers (each 2 BLAS threads, contiguous chromosome-sorted chunks), combine.
set -euo pipefail
REPO=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
OUT=$REPO/Analysis/Hypotheses/H27-20261008T230909Z
SC=$REPO/scripts/analysis/hypotheses_2026q4/H27_kidney_hf_coloc
RS=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
FTP=https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics
NW=4

export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2
export MASLD_PROJECT_ROOT=$REPO
unset LD_PANEL UKBB_LD_DIR || true   # release default EUR panel (PolyFun UKB)

echo "host $(hostname)  start $(date -u +%FT%TZ)  job ${SLURM_JOB_ID:-NA}"
test -f "$OUT/PRESPEC_H27.md"
test ! -e "$OUT/results"   # never overwrite an earlier run
sha256sum "$OUT/PRESPEC_H27.md" "$SC/h27_kidney_hf_coloc.R" "$SC/run_h27.sbatch" \
  | tee "$OUT/logs/script_sha256_${SLURM_JOB_ID:-NA}.txt"

# ---- outcome downloads (GWAS Catalog FTP; each file in its own directory) ----
fetch() {  # acc subdir file has_md5
  local acc=$1 sub=$2 file=$3 has_md5=$4 dir=$OUT/inputs/outcomes/$1
  mkdir -p "$dir"
  if [[ -f "$dir/sha256.txt" ]] && (cd "$dir" && sha256sum -c --quiet sha256.txt); then
    echo "$acc present, sha256 ok"; return 0
  fi
  curl -fsS --retry 5 --retry-delay 30 -o "$dir/$file" "$FTP/$sub/$acc/$file"
  if [[ $has_md5 == 1 ]]; then
    curl -fsS --retry 5 -o "$dir/md5sum.txt" "$FTP/$sub/$acc/md5sum.txt"
    local want got
    want=$(awk -v f="$file" '$2==f || $2=="*"f {print $1}' "$dir/md5sum.txt")
    got=$(md5sum "$dir/$file" | cut -d' ' -f1)
    echo "$acc md5 published $want computed $got"
    [[ -n $want && $want == "$got" ]]
  else
    echo "$acc: no md5 published; sha256 recorded only"
  fi
  (cd "$dir" && sha256sum "$file" > sha256.txt && cat sha256.txt)
}
fetch GCST90103634 GCST90103001-GCST90104000 GCST90103634_buildGRCh37.tsv 1
fetch GCST008065   GCST008001-GCST009000     CKD_overall_EA_JW_20180223_nstud23.dbgap.txt.gz 0
fetch GCST90728695 GCST90728001-GCST90729000 GCST90728695.tsv.gz 1

# ---- prep ----
"$RS" --vanilla "$SC/h27_kidney_hf_coloc.R" prep "$REPO" "$OUT"

# ---- coloc workers ----
pids=()
for ((w = 0; w < NW; w++)); do
  "$RS" --vanilla "$SC/h27_kidney_hf_coloc.R" coloc "$REPO" "$OUT" "$w" "$NW" \
    > "$OUT/logs/coloc_worker${w}_${SLURM_JOB_ID:-NA}.log" 2>&1 &
  pids+=($!)
done
fail=0
for p in "${pids[@]}"; do wait "$p" || fail=1; done
if [[ $fail == 1 ]]; then echo "a coloc worker failed; see logs"; exit 1; fi

# ---- combine (calls, published-list overlap, controls, permutation null) ----
"$RS" --vanilla "$SC/h27_kidney_hf_coloc.R" combine "$REPO" "$OUT"
echo "end $(date -u +%FT%TZ)"
sstat -j "${SLURM_JOB_ID:-0}.batch" --format=MaxRSS 2>/dev/null || true
