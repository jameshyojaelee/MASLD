#!/bin/bash
#SBATCH --job-name=pazoki_derived_fix
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/orientation_checks/logs/slurm-%j.out

# Step 5 of the Pazoki 2022 PDFF (GCST90267352) orientation fix. Step 3 corrected the four primary copies; this job
# corrects the copies derived from them:
#   sumstats_af/2022_36402844_PDFF_EUR_af.tsv                    effect/other labels exchanged
#   lead_snps/..._preprocessed_leadSNPs_merged.tsv (11 rows)     effect/other labels exchanged
#   sumstats/Pazoki_PDFF_reformatted_hg19.tsv                    labels exchanged, position restored
#   MR_Data/hg19/GCST90267352_PDFF_Pazoki2022_hg19.tsv.gz        hg19_chr, hg19_pos restored
#   config/gwas_af_sources.tsv                                   sha256 of the _af.tsv copy
# The source positions are GRCh37; the old hg19 columns lifted them a second time. beta, se, p and every other
# column stay byte-identical. Corrected copies are written and validated first. Only if every check passes is each
# canonical file replaced by one atomic rename, with the previous file kept beside it under a hard link.
set -euo pipefail
trap 'echo "FAILED at line $LINENO (stage ${STAGE:-start})"' ERR

REPO=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
G=$REPO/GWAS
SC=$REPO/scripts/analysis/gwas_orientation
PY=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python
FIX="$PY -I $SC/fix_pazoki_derived_copies.py"
FIRST=$G/orientation_checks/pazoki_pdff-20261010T022251Z
UTC=$(date -u +%Y%m%dT%H%M%SZ)
OUT=$G/orientation_checks/pazoki_pdff_derived_copies-$UTC
echo "host $(hostname)  start $(date -u +%FT%TZ)  job ${SLURM_JOB_ID:-NA}  out $OUT"

A=$G/MR_Data/GCST90267352_PDFF_Pazoki2022.tsv.gz
MRHG19=$G/MR_Data/hg19/GCST90267352_PDFF_Pazoki2022_hg19.tsv.gz
C=$G/finemapping/data/sumstats/2022_36402844_PDFF_EUR_preprocessed.tsv
R1=$G/finemapping/data/sumstats/2021_34128465_PDFF_EUR_preprocessed.tsv
R2=$G/finemapping/data/sumstats/2021_34957434_PDFF_EUR_preprocessed.tsv
AF=$G/finemapping/data/sumstats_af/2022_36402844_PDFF_EUR_af.tsv
LEAD=$G/finemapping/data/lead_snps/2022_36402844_PDFF_EUR_preprocessed_leadSNPs_merged.tsv
REFHG19=$G/finemapping/data/sumstats/Pazoki_PDFF_reformatted_hg19.tsv
SOURCES=$G/finemapping/config/gwas_af_sources.tsv
STUDY=2022_36402844_PDFF_EUR

# values recorded by step 3 (swap_log.tsv) and in config/gwas_af_sources.tsv before this job
A_MD5=2d56a713f4bd2b4b6e2a94278482bdf5
B_MD5=028cdb05dc1a6aa12b2b2ea6c6f0ef96
C_MD5=1389da84f9607dc49757dff5c3a7e0f2
AF_SHA_OLD=3f0942c9d1cb7cf99f01ac2e39c634be1ff70f84225ef3d0d2aff64f818990a7
B_ROWS=8975721          # hg19_summary.json from step 1
B_POS_CHANGED=8895197

md5_of() { md5sum "$1" | cut -d' ' -f1; }
STAGE=preconditions
[[ $(md5_of "$A") == "$A_MD5" ]] || { echo "corrected MR_Data copy changed"; exit 2; }
[[ $(md5_of "$MRHG19") == "$B_MD5" ]] || { echo "MR_Data hg19 copy changed since step 3"; exit 2; }
[[ $(md5_of "$C") == "$C_MD5" ]] || { echo "corrected preprocessed copy changed"; exit 2; }
[[ $(sha256sum "$AF" | cut -d' ' -f1) == "$AF_SHA_OLD" ]] || { echo "_af.tsv is not the file hashed in gwas_af_sources.tsv"; exit 2; }
[[ $(awk -F'\t' -v s="$STUDY" 'NR>1 && $1==s {print $11}' "$SOURCES") == "$AF_SHA_OLD" ]] || { echo "sources line does not hold the expected hash"; exit 2; }
test ! -e "$OUT"
mkdir "$OUT" "$OUT/corrected" "$OUT/work"
sha256sum "$SC/fix_pazoki_derived_copies.py" "$SC/run_fix_derived_copies.sbatch" | tee "$OUT/script_sha256.txt"
$PY -I --version > "$OUT/environment.txt" 2>&1
$PY -I -m pip freeze > "$OUT/pip_freeze.txt"
{ printf '%s  AF  %s\n' "$(md5_of "$AF")" "$AF"
  printf '%s  LEAD  %s\n' "$(md5_of "$LEAD")" "$LEAD"
  printf '%s  REFHG19  %s\n' "$(md5_of "$REFHG19")" "$REFHG19"
  printf '%s  MRHG19  %s\n' "$(md5_of "$MRHG19")" "$MRHG19"
  printf '%s  SOURCES  %s\n' "$(md5_of "$SOURCES")" "$SOURCES"; } > "$OUT/md5_before.txt"
cat "$OUT/md5_before.txt"
STAGE=selftest
$FIX selftest

reader() { if [[ $1 == *.gz ]]; then echo zcat; else echo cat; fi; }
colidx() ( set +o pipefail; $(reader "$1") "$1" | head -n1 | tr '\t' '\n' | grep -n -x -F "$2" | cut -d: -f1 )
swap_awk() {  # c1 c2: exchange two tab-separated columns, leave every other byte alone
  awk -F'\t' -v OFS='\t' -v c1="$1" -v c2="$2" '
    NR==1 {nf=NF; print; next}
    NF!=nf {printf "field count %d != %d at line %d\n", NF, nf, NR > "/dev/stderr"; exit 1}
    {t=$c1; $c1=$c2; $c2=t; print}'
}
VLOG=$OUT/shell_validation.tsv
printf 'file\tcheck\tok\n' > "$VLOG"
ALL_OK=1
record() { printf '%s\t%s\t%s\n' "$1" "$2" "$3" >> "$VLOG"; [[ $3 == 1 ]] || ALL_OK=0; }
flag() { [[ $1 == "$2" ]] && echo 1 || echo 0; }

exchange_alleles() {  # label in out n1 n2: write out with the two columns exchanged and validate it against in
  local label=$1 in=$2 out=$3 n1=$4 n2=$5 c1 c2
  c1=$(colidx "$in" "$n1"); c2=$(colidx "$in" "$n2")
  [[ $c1 =~ ^[0-9]+$ && $c2 =~ ^[0-9]+$ ]] || { echo "cannot find $n1/$n2 in $in"; exit 4; }
  echo "$label: exchanging columns $c1 ($n1) and $c2 ($n2)"
  cat "$in" | swap_awk "$c1" "$c2" > "$out"
  chmod --reference="$in" "$out"
  local l_in l_out o_in o_out rt_in rt_out sw1 sw2
  l_in=$(wc -l < "$in"); l_out=$(wc -l < "$out")
  o_in=$(cut --complement -f"$c1,$c2" "$in" | md5sum | cut -d' ' -f1)
  o_out=$(cut --complement -f"$c1,$c2" "$out" | md5sum | cut -d' ' -f1)
  rt_in=$(md5sum < "$in" | cut -d' ' -f1)
  rt_out=$(swap_awk "$c1" "$c2" < "$out" | md5sum | cut -d' ' -f1)
  sw1=$(paste <(cut -f"$c1" "$in") <(cut -f"$c2" "$out") | awk -F'\t' 'NR>1 && $1!=$2{bad++} END{print bad+0}')
  sw2=$(paste <(cut -f"$c2" "$in") <(cut -f"$c1" "$out") | awk -F'\t' 'NR>1 && $1!=$2{bad++} END{print bad+0}')
  record "$label" "row count equal" "$(flag "$l_in" "$l_out")"
  record "$label" "all columns other than $n1/$n2 byte-identical" "$(flag "$o_in" "$o_out")"
  record "$label" "swapping back reproduces the original text exactly" "$(flag "$rt_in" "$rt_out")"
  record "$label" "corrected $n2 column equals original $n1 column on every row" "$(flag "$sw1" 0)"
  record "$label" "corrected $n1 column equals original $n2 column on every row" "$(flag "$sw2" 0)"
  echo "$label: lines $l_in -> $l_out  $(date -u +%FT%TZ)"
}

# 1. _af.tsv
STAGE=af
exchange_alleles af "$AF" "$OUT/corrected/$(basename "$AF")" allele1 allele2
# informational: the _af copy should now equal the corrected preprocessed copy in its first seven columns
{ printf 'check\tvalue\n'
  printf 'af cut -f1-7 md5 equals corrected C file md5\t%s\n' "$(flag "$(cut -f1-7 "$OUT/corrected/$(basename "$AF")" | md5sum | cut -d' ' -f1)" "$C_MD5")"
  printf 'af allele columns differing from corrected C, row by row\t%s\n' "$(paste <(cut -f3,4 "$OUT/corrected/$(basename "$AF")") <(cut -f3,4 "$C") | awk -F'\t' 'NR>1 && ($1!=$3||$2!=$4){bad++} END{print bad+0}')"
} > "$OUT/af_vs_corrected_C.tsv"
cat "$OUT/af_vs_corrected_C.tsv"

# 2. lead-SNP table
STAGE=lead
exchange_alleles lead "$LEAD" "$OUT/corrected/$(basename "$LEAD")" allele1 allele2
$FIX lead-check "$LEAD" "$OUT/corrected/$(basename "$LEAD")" "$C" "$OUT/lead_check.tsv"

# 3. reformatted hg19 copy: labels, then position
STAGE=reformatted
exchange_alleles reformatted_labels "$REFHG19" "$OUT/work/E_swapped.tsv" allele1 allele2
EOUT=$OUT/corrected/$(basename "$REFHG19")
$FIX positions "$MRHG19" "$OUT/work/E_swapped.tsv" "$EOUT" "$OUT/positions_report.json" "$FIRST/hg19_offsets_by_chr.tsv"
chmod --reference="$REFHG19" "$EOUT"
record reformatted "row count equal" "$(flag "$(wc -l < "$REFHG19")" "$(wc -l < "$EOUT")")"
record reformatted "columns chromosome, beta, se, pval byte-identical" \
  "$(flag "$(cut --complement -f2,3,4 "$REFHG19" | md5sum)" "$(cut --complement -f2,3,4 "$EOUT" | md5sum)")"
record reformatted "final allele2 equals original allele1 on every row" \
  "$(flag "$(paste <(cut -f3 "$REFHG19") <(cut -f4 "$EOUT") | awk -F'\t' 'NR>1 && $1!=$2{bad++} END{print bad+0}')" 0)"
record reformatted "final allele1 equals original allele2 on every row" \
  "$(flag "$(paste <(cut -f4 "$REFHG19") <(cut -f3 "$EOUT") | awk -F'\t' 'NR>1 && $1!=$2{bad++} END{print bad+0}')" 0)"
read -r PR_OK <<<"$($PY -I -c "import json,sys; d=json.load(open(sys.argv[1])); print(int(all(v == 0 for v in d['pairing_mismatches'].values()) and d['position_changed'] + d['position_unchanged'] == d['rows_E']))" "$OUT/positions_report.json")"
record reformatted "every row paired with the MR_Data hg19 row on chromosome, alleles, beta, se, p and old position" "$PR_OK"

# 4. MR_Data hg19 copy: hg19_chr and hg19_pos
STAGE=mr_hg19
BOUT=$OUT/corrected/$(basename "$MRHG19")
$FIX hg19-columns "$MRHG19" "$BOUT" "$OUT/hg19_columns_report.json"
chmod --reference="$MRHG19" "$BOUT"
record mr_hg19 "row count equal" "$(flag "$(zcat "$MRHG19" | wc -l)" "$(zcat "$BOUT" | wc -l)")"
record mr_hg19 "columns 1 to 10 byte-identical" \
  "$(flag "$(zcat "$MRHG19" | cut -f1-10 | md5sum)" "$(zcat "$BOUT" | cut -f1-10 | md5sum)")"
record mr_hg19 "hg19_chr equals chromosome and hg19_pos equals base_pair_location on every row" \
  "$(flag "$(zcat "$BOUT" | awk -F'\t' 'NR>1 && ($11!=$3||$12!=$4){bad++} END{print bad+0}')" 0)"
read -r HG_ROWS HG_POS <<<"$($PY -I -c "import json,sys; d=json.load(open(sys.argv[1])); print(d['rows'], d['hg19_pos_changed'])" "$OUT/hg19_columns_report.json")"
record mr_hg19 "rows and changed hg19_pos equal the step 1 counts ($B_ROWS, $B_POS_CHANGED)" \
  "$(flag "$HG_ROWS $HG_POS" "$B_ROWS $B_POS_CHANGED")"

record reformatted "position equals hg19_pos of the corrected MR_Data hg19 copy on every row" \
  "$(flag "$(paste <(cut -f2 "$EOUT") <(zcat "$BOUT" | cut -f12) | awk -F'\t' 'NR>1 && $1!=$2{bad++} END{print bad+0}')" 0)"

# 5. positions against the two other PDFF studies (GRCh37)
STAGE=agreement
$FIX agreement "$OUT/agreement.json" "$A" "$MRHG19" "$BOUT" "$REFHG19" "$EOUT" "$R1" "$R2"

cat "$VLOG"
STAGE=gate
if [[ $ALL_OK == 1 ]]; then : > "$OUT/VALIDATION_PASS"; echo VALIDATION_PASS; else echo VALIDATION_FAIL; exit 5; fi

# 6. replace the canonical files (one atomic rename each; the previous file stays beside it as a hard link)
STAGE=replace
printf 'key\tcanonical\tbackup\tmd5_before\tmd5_backup\tmd5_corrected\tmd5_canonical_after\tperm_before\tperm_after\n' > "$OUT/swap_log.tsv"
replace() {  # key canonical corrected_copy backup_tag
  local k=$1 canon=$2 src=$3 tag=$4 backup tmp want now corr perm_before after
  backup="$canon.${tag}_$UTC"; tmp="$canon.new_$UTC"
  want=$(awk -v k="$k" '$2==k{print $1}' "$OUT/md5_before.txt")
  test -f "$src"; test ! -e "$backup"; test ! -e "$tmp"
  now=$(md5_of "$canon")
  [[ $now == "$want" ]] || { echo "$k: canonical file changed since validation ($now != $want); not replacing"; exit 3; }
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
  [[ $after == "$corr" ]] || { echo "$k: canonical md5 after the rename differs from the corrected file"; exit 6; }
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$k" "$canon" "$backup" "$want" "$(md5_of "$backup")" "$corr" "$after" \
    "$perm_before" "$(stat -c '%U:%G %a' "$canon")" >> "$OUT/swap_log.tsv"
  echo "$k replaced"
}
replace AF "$AF" "$OUT/corrected/$(basename "$AF")" before_label_fix
replace LEAD "$LEAD" "$OUT/corrected/$(basename "$LEAD")" before_label_fix
replace REFHG19 "$REFHG19" "$EOUT" before_label_and_position_fix
replace MRHG19 "$MRHG19" "$BOUT" before_hg19_position_fix

# 7. the hash of the _af.tsv copy in the AF source list (one field of one line)
STAGE=sources
AF_SHA_NEW=$(sha256sum "$AF" | cut -d' ' -f1)
$FIX update-sources "$SOURCES" "$SOURCES.new_$UTC" "$STUDY" "$AF_SHA_OLD" "$AF_SHA_NEW"
[[ $(diff "$SOURCES" "$SOURCES.new_$UTC" | grep -c '^<') == 1 ]] || { echo "sources diff is not one line"; exit 7; }
cp -p "$SOURCES" "$SOURCES.before_af_hash_update_$UTC"
chmod --reference="$SOURCES" "$SOURCES.new_$UTC"
mv -f "$SOURCES.new_$UTC" "$SOURCES"
diff "$SOURCES.before_af_hash_update_$UTC" "$SOURCES" | cut -c1-200 | tee "$OUT/sources_diff.txt" || true

# 8. record
STAGE=record
{ for k in AF LEAD REFHG19 MRHG19 SOURCES; do
    p=$(awk -v k="$k" '$2==k{print $3}' "$OUT/md5_before.txt")
    printf '%s  %s  %s\n' "$(sha256sum "$p" | cut -d' ' -f1)" "$k" "$p"
  done; } > "$OUT/sha256_after.txt"
: > "$OUT/REPLACED_OK"
cat "$OUT/swap_log.tsv" "$OUT/sha256_after.txt"
echo "end $(date -u +%FT%TZ)"
