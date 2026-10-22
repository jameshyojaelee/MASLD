#!/bin/bash
#SBATCH --job-name=pazoki_orient_fix
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=96G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/orientation_checks/logs/slurm-%j.out

# Step 2: write corrected copies of the four Pazoki PDFF files (effect/other labels exchanged; beta, z, se, p, n and
# every other column untouched) into $ORIENT_DIR/corrected and validate them. No input file is modified.
# Requires ORIENT_DIR (the step 1 results directory). Stops before writing anything if the re-check of the
# as-downloaded files does not choose the label swap.
set -euo pipefail
: "${ORIENT_DIR:?ORIENT_DIR is required}"
REPO=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
G=$REPO/GWAS
SC=$REPO/scripts/analysis/gwas_orientation
PY=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python
CD=$ORIENT_DIR/corrected
echo "host $(hostname)  start $(date -u +%FT%TZ)  job ${SLURM_JOB_ID:-NA}  orient_dir $ORIENT_DIR"
test -f "$ORIENT_DIR/decision.json"
test ! -e "$CD"
test ! -e "$ORIENT_DIR/recheck_original"
test ! -e "$ORIENT_DIR/corrected_validation"
sha256sum "$SC/pazoki_orientation_check.py" "$SC/run_make_corrected.sbatch" | tee "$ORIENT_DIR/script_sha256_step2.txt"

declare -A SRC=(
  [GCST90267352_PDFF_Pazoki2022.tsv.gz]="$G/MR_Data/GCST90267352_PDFF_Pazoki2022.tsv.gz"
  [GCST90267352_PDFF_Pazoki2022_hg19.tsv.gz]="$G/MR_Data/hg19/GCST90267352_PDFF_Pazoki2022_hg19.tsv.gz"
  [2022_36402844_PDFF_EUR_preprocessed.tsv]="$G/finemapping/data/sumstats/2022_36402844_PDFF_EUR_preprocessed.tsv"
  [pdff.sumstats.gz]="$G/ldsc/munged/pdff.sumstats.gz")
declare -A NAMES=(
  [GCST90267352_PDFF_Pazoki2022.tsv.gz]="effect_allele other_allele"
  [GCST90267352_PDFF_Pazoki2022_hg19.tsv.gz]="effect_allele other_allele"
  [2022_36402844_PDFF_EUR_preprocessed.tsv]="allele1 allele2"
  [pdff.sumstats.gz]="A1 A2")

reader() { if [[ $1 == *.gz ]]; then echo zcat; else echo cat; fi; }
colidx() ( set +o pipefail; $(reader "$1") "$1" | head -n1 | tr '\t' '\n' | grep -n -x -F "$2" | cut -d: -f1 )
swap_awk() {  # c1 c2: exchange two tab-separated columns, leave every other byte alone
  awk -F'\t' -v OFS='\t' -v c1="$1" -v c2="$2" '
    NR==1 {nf=NF; print; next}
    NF!=nf {printf "field count %d != %d at line %d\n", NF, nf, NR > "/dev/stderr"; exit 1}
    {t=$c1; $c1=$c2; $c2=t; print}'
}

# 1. re-check the as-downloaded files with the key-based derived-file comparisons
RC=$ORIENT_DIR/recheck_original
mkdir "$RC"
"$PY" -I "$SC/pazoki_orientation_check.py" "$RC"
grep -q '"fix": "swap_effect_other_labels_keep_beta_z_se_p"' "$RC/decision.json" || { echo "STOP: re-check did not choose the label swap"; exit 3; }

# 2. write the corrected copies
mkdir "$CD"
: > "$ORIENT_DIR/shell_validation.tsv"; printf 'file\tcheck\tok\n' >> "$ORIENT_DIR/shell_validation.tsv"
ALL_OK=1
for f in "${!SRC[@]}"; do
  in=${SRC[$f]}; out=$CD/$f
  read -r n1 n2 <<<"${NAMES[$f]}"
  c1=$(colidx "$in" "$n1"); c2=$(colidx "$in" "$n2")
  [[ -n $c1 && -n $c2 && $c1 =~ ^[0-9]+$ && $c2 =~ ^[0-9]+$ ]] || { echo "cannot find $n1/$n2 in $in"; exit 4; }
  echo "$f: exchanging columns $c1 ($n1) and $c2 ($n2)"
  rd=$(reader "$in")
  if [[ $in == *.gz ]]; then
    $rd "$in" | swap_awk "$c1" "$c2" | gzip -n -6 > "$out"
  else
    $rd "$in" | swap_awk "$c1" "$c2" > "$out"
  fi
  chmod --reference="$in" "$out"
  # 3. byte-level validation
  l_in=$($rd "$in" | wc -l); l_out=$($rd "$out" | wc -l)
  o_in=$($rd "$in" | cut --complement -f"$c1,$c2" | md5sum | cut -d' ' -f1)
  o_out=$($rd "$out" | cut --complement -f"$c1,$c2" | md5sum | cut -d' ' -f1)
  rt_in=$($rd "$in" | md5sum | cut -d' ' -f1)
  rt_out=$($rd "$out" | swap_awk "$c1" "$c2" | md5sum | cut -d' ' -f1)
  sw=$(paste <($rd "$in" | cut -f"$c1") <($rd "$out" | cut -f"$c2") | awk -F'\t' 'NR>1 && $1!=$2{bad++} END{print bad+0}')
  for chk in "row count equal:$([[ $l_in == "$l_out" ]] && echo 1 || echo 0)" \
             "all columns other than $n1/$n2 byte-identical:$([[ $o_in == "$o_out" ]] && echo 1 || echo 0)" \
             "swapping back reproduces the original text exactly:$([[ $rt_in == "$rt_out" ]] && echo 1 || echo 0)" \
             "corrected $n2 column equals original $n1 column on every row:$([[ $sw == 0 ]] && echo 1 || echo 0)"; do
    ok=${chk##*:}; name=${chk%:*}
    printf '%s\t%s\t%s\n' "$f" "$name" "$ok" >> "$ORIENT_DIR/shell_validation.tsv"
    [[ $ok == 1 ]] || ALL_OK=0
  done
  echo "$f: lines $l_in -> $l_out; written $(date -u +%FT%TZ)"
done
cat "$ORIENT_DIR/shell_validation.tsv"

# 4. rerun the full orientation check on the corrected copies
VAL=$ORIENT_DIR/corrected_validation
mkdir "$VAL"
"$PY" -I "$SC/pazoki_orientation_check.py" "$VAL" corrected "$CD"
grep -q '"corrected_ok": true' "$VAL/decision.json" || ALL_OK=0
if [[ $ALL_OK == 1 ]]; then
  : > "$ORIENT_DIR/VALIDATION_PASS"; echo VALIDATION_PASS
else
  echo VALIDATION_FAIL; exit 5
fi
echo "end $(date -u +%FT%TZ)"
