#!/bin/bash
#SBATCH --job-name=pazoki_orient_swap
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/orientation_checks/logs/slurm-%j.out

# Step 3: put the validated corrected copies at the four canonical paths and keep each original beside it as
# <name>.as_downloaded_<UTC>. Each canonical file is replaced by one atomic rename; the original survives under
# its backup name (hard link made first), so no path is ever missing and nothing is deleted.
# Requires ORIENT_DIR with VALIDATION_PASS (written by run_make_corrected.sbatch) and md5_before.txt.
set -euo pipefail
: "${ORIENT_DIR:?ORIENT_DIR is required}"
G=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS
CD=$ORIENT_DIR/corrected
test -f "$ORIENT_DIR/VALIDATION_PASS"
test -f "$ORIENT_DIR/md5_before.txt"
test ! -e "$ORIENT_DIR/swap_log.tsv"
UTC=$(date -u +%Y%m%dT%H%M%SZ)
echo "host $(hostname)  start $(date -u +%FT%TZ)  job ${SLURM_JOB_ID:-NA}  utc $UTC"
sha256sum /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/gwas_orientation/run_swap_in_place.sbatch | tee "$ORIENT_DIR/script_sha256_step3.txt"

declare -A CANON=(
  [A]="$G/MR_Data/GCST90267352_PDFF_Pazoki2022.tsv.gz"
  [B]="$G/MR_Data/hg19/GCST90267352_PDFF_Pazoki2022_hg19.tsv.gz"
  [C]="$G/finemapping/data/sumstats/2022_36402844_PDFF_EUR_preprocessed.tsv"
  [D]="$G/ldsc/munged/pdff.sumstats.gz")
md5_of() { md5sum "$1" | cut -d' ' -f1; }
printf 'key\tcanonical\tbackup\tmd5_original\tmd5_backup\tmd5_corrected\tmd5_canonical_after\tperm_before\tperm_after\n' > "$ORIENT_DIR/swap_log.tsv"
: > "$ORIENT_DIR/md5_after.txt"
: > "$ORIENT_DIR/sha256_after.txt"   # sha256 of each original (backup) and each corrected canonical file, to look up in sha256 manifests
for k in A B C D; do
  canon=${CANON[$k]}; src=$CD/$(basename "$canon"); backup="$canon.as_downloaded_$UTC"; tmp="$canon.new_$UTC"
  want=$(awk -v k="$k" '$2==k{print $1}' "$ORIENT_DIR/md5_before.txt")
  [[ -n $want ]] || { echo "no recorded md5 for $k"; exit 2; }
  test -f "$src"
  test ! -e "$backup"
  test ! -e "$tmp"
  now=$(md5_of "$canon")
  [[ $now == "$want" ]] || { echo "$k: canonical file changed since verification ($now != $want); not replacing"; exit 3; }
  corr=$(md5_of "$src")
  perm_before=$(stat -c '%U:%G %a' "$canon")
  cp "$src" "$tmp"
  chmod --reference="$canon" "$tmp"
  chgrp --reference="$canon" "$tmp" 2>/dev/null || true
  [[ $(md5_of "$tmp") == "$corr" ]] || { echo "$k: copy differs from corrected file"; exit 4; }
  ln "$canon" "$backup"
  [[ $(md5_of "$backup") == "$want" ]] || { echo "$k: backup md5 mismatch"; exit 5; }
  mv -f "$tmp" "$canon"
  after=$(md5_of "$canon")
  [[ $after == "$corr" ]] || { echo "$k: canonical md5 after swap differs from corrected"; exit 6; }
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$k" "$canon" "$backup" "$want" "$(md5_of "$backup")" "$corr" "$after" "$perm_before" "$(stat -c '%U:%G %a' "$canon")" >> "$ORIENT_DIR/swap_log.tsv"
  printf '%s  %s  %s\n%s  %s  %s\n' "$want" "$k.original" "$backup" "$after" "$k.corrected" "$canon" >> "$ORIENT_DIR/md5_after.txt"
  printf '%s  %s.original  %s\n%s  %s.corrected  %s\n' "$(sha256sum "$backup" | cut -d' ' -f1)" "$k" "$backup" "$(sha256sum "$canon" | cut -d' ' -f1)" "$k" "$canon" >> "$ORIENT_DIR/sha256_after.txt"
  echo "$k done"
done
cat "$ORIENT_DIR/swap_log.tsv"
echo "end $(date -u +%FT%TZ)"
