#!/bin/bash
#SBATCH --job-name=sldsc_download
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=io
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=GWAS/ldsc/program_heritability_20260830/logs/10_download_%j.out
#SBATCH --error=GWAS/ldsc/program_heritability_20260830/logs/10_download_%j.err
#
# Fetch the three reference files that stratified LD score regression needs and
# that this project does not already have. Everything else is already on disk:
# the 1000 Genomes Phase 3 European PLINK panel and the HapMap3 SNP lists ship
# with the gs-MAP resource under data/gsmap_resource/.
#
# The Broad's old alkesgroup web directory no longer serves these files. The
# files are now distributed from Zenodo record 10515792, "S-LDSC reference
# files", DOI 10.5281/zenodo.10515792, deposited 2024-01-15 under CC-BY-4.0.
# The MD5 sums below were read from the Zenodo record on 2026-08-30 and are
# checked after download; a mismatch stops the job.
#
# Total download: 775 MB compressed. Roughly 6 GB after extraction.

set -euo pipefail
cd "${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"

DEST=GWAS/ldsc_ref/sldsc

# Refuse to run if another copy of this job is already downloading. On
# 2026-08-30 two of these ran at once on different nodes and both wrote the same
# tarball with "curl -o", which truncates. The files survived that time, because
# both were fetching identical bytes from the same URL, but the extraction that
# followed was interrupted and left a tree whose file counts looked complete.
# One lock removes the whole class of problem.
mkdir -p "$(dirname "$DEST")"
exec 9>"$(dirname "$DEST")/.sldsc_download.lock"
if ! flock -n 9; then
  echo "another sldsc_download job holds the lock; refusing to run concurrently" >&2
  exit 1
fi

mkdir -p "$DEST/tarballs"

ZEN=https://zenodo.org/api/records/10515792/files

fetch () {  # fetch <filename> <md5>
  local name=$1 md5=$2 out="$DEST/tarballs/$1"
  if [[ -f $out ]] && [[ "$(md5sum "$out" | cut -d' ' -f1)" == "$md5" ]]; then
    echo "already have $name"
    return
  fi
  echo ">> downloading $name"
  curl -fL --retry 5 --retry-delay 20 -o "$out" "$ZEN/$name/content"
  local got
  got=$(md5sum "$out" | cut -d' ' -f1)
  if [[ "$got" != "$md5" ]]; then
    echo "MD5 mismatch for $name: expected $md5, got $got" >&2
    exit 1
  fi
  echo "   md5 ok"
}

fetch 1000G_Phase3_baselineLD_v2.2_ldscores.tgz b261e0caf06a003e7522938e01b3d349
fetch 1000G_Phase3_frq.tgz                      ac29686ffd5b6378789857a522ebca77
fetch 1000G_Phase3_weights_hm3_no_MHC.tgz       a98ac0f089ee285177544a3e6e721ca3

# Extract into a staging directory and move it into place only once every
# archive has been unpacked. An interrupted tar otherwise leaves files that are
# present but truncated, and a count of 22 per chromosome looks complete when it
# is not.
already_extracted () {
  local d n
  for d in baselineLD_v2.2 1000G_Phase3_frq 1000G_Phase3_weights_hm3_no_MHC; do
    [[ -d "$DEST/$d" ]] || return 1
  done
  n=$(find "$DEST/baselineLD_v2.2" -name 'baselineLD.*.l2.ldscore.gz' | wc -l)
  [[ "$n" -eq 22 ]] || return 1
  n=$(find "$DEST" -name '1000G.EUR.QC.*.frq' | wc -l)
  [[ "$n" -eq 22 ]] || return 1
  n=$(find "$DEST" -name 'weights.hm3_noMHC.*.l2.ldscore.gz' | wc -l)
  [[ "$n" -eq 22 ]] || return 1
  return 0
}

if already_extracted; then
  echo ">> reference already unpacked with 22 chromosomes of each; verifying instead of re-extracting"
  bad=0
  while IFS= read -r g; do
    gzip -t "$g" 2>/dev/null || { echo "corrupt: $g"; bad=1; }
  done < <(find "$DEST/baselineLD_v2.2" "$DEST/1000G_Phase3_weights_hm3_no_MHC" -name "*.gz")
  [[ $bad -eq 1 ]] && exit 1
  echo "   every extracted archive passes gzip -t"
else

echo ">> extracting into staging"
# Nothing here deletes recursively. Anything in the way is moved aside with a
# timestamp so it stays recoverable, and the user decides when to reclaim it.
STAMP=$(date +%Y%m%dT%H%M%S)
STAGE="$DEST/.staging"
[[ -e $STAGE ]] && mv "$STAGE" "$DEST/.staging.superseded-$STAMP"
mkdir -p "$STAGE/baselineLD_v2.2"
tar -xzf "$DEST/tarballs/1000G_Phase3_baselineLD_v2.2_ldscores.tgz" -C "$STAGE/baselineLD_v2.2" --strip-components=0
tar -xzf "$DEST/tarballs/1000G_Phase3_frq.tgz"                      -C "$STAGE"
tar -xzf "$DEST/tarballs/1000G_Phase3_weights_hm3_no_MHC.tgz"       -C "$STAGE"

echo ">> testing every extracted gzip before moving into place"
bad=0
while IFS= read -r g; do
  gzip -t "$g" 2>/dev/null || { echo "corrupt after extraction: $g"; bad=1; }
done < <(find "$STAGE" -name "*.gz")
[[ $bad -eq 1 ]] && exit 1

echo ">> moving staging into place"
for d in baselineLD_v2.2 1000G_Phase3_frq 1000G_Phase3_weights_hm3_no_MHC; do
  if [[ -d "$STAGE/$d" ]]; then
    [[ -e "$DEST/$d" ]] && mv "$DEST/$d" "$DEST/$d.superseded-$STAMP"
    mv "$STAGE/$d" "$DEST/$d"
  fi
done
rmdir "$STAGE" 2>/dev/null || true
fi

echo ">> checking that every chromosome arrived"
missing=0
for c in $(seq 1 22); do
  for f in \
    "$(find "$DEST/baselineLD_v2.2" -name "baselineLD.$c.l2.ldscore.gz" | head -1)" \
    "$(find "$DEST/baselineLD_v2.2" -name "baselineLD.$c.annot.gz"      | head -1)" \
    "$(find "$DEST/baselineLD_v2.2" -name "baselineLD.$c.l2.M_5_50"     | head -1)" \
    "$(find "$DEST" -name "1000G.EUR.QC.$c.frq"                          | head -1)" \
    "$(find "$DEST" -name "weights.hm3_noMHC.$c.l2.ldscore.gz"           | head -1)" ; do
    [[ -z $f ]] && { echo "MISSING a file for chromosome $c"; missing=1; }
  done
done
[[ $missing -eq 1 ]] && exit 1

echo ">> baselineLD annotation count"
# Count the header fields without ever closing the pipe early. Anything that
# stops reading after the first line - head -1, or awk with NR>1{exit} - makes
# zcat die on SIGPIPE and return 141, which under `set -o pipefail` fails the
# job even though the correct answer was already printed. Reading the file to
# EOF costs a few seconds and removes the failure mode entirely.
LD22=$(find "$DEST/baselineLD_v2.2" -name 'baselineLD.22.l2.ldscore.gz' -print -quit)
n_baseline_annot=$(zcat "$LD22" | awk 'NR==1{print NF-3; f=1} END{if(!f) print 0}')
echo "$n_baseline_annot annotations in the baselineLD v2.2 model"
if [[ "$n_baseline_annot" -lt 90 ]]; then
  echo "baselineLD v2.2 should carry 97 annotations; found $n_baseline_annot" >&2
  exit 1
fi

echo ">> checking that baselineLD sits on the same SNPs as our PLINK panel"
# The two files must hold the SAME chr22 SNPs, not merely overlapping ones, so
# the test is set identity on the SNP identifiers and it fails on any
# difference. A coverage floor would have been the wrong test: at 50 percent it
# passes a reference that has silently lost half of chromosome 22, which is
# exactly the failure it exists to catch.
#
# Note this applies to baselineLD.<c>.annot.gz, which carries every reference
# SNP. The sibling baselineLD.<c>.l2.ldscore.gz is a different and much smaller
# set, the HapMap3 regression SNPs, and there equality would be wrong: LDSC
# intersects those with the summary statistics at regression time.
PLINK=data/gsmap_resource/LD_Reference_Panel/1000G_EUR_Phase3_plink
ANNOT22=$(find "$DEST/baselineLD_v2.2" -name 'baselineLD.22.annot.gz' -print -quit)
if [[ -z $ANNOT22 ]]; then
  echo "baselineLD.22.annot.gz not found under $DEST/baselineLD_v2.2" >&2
  exit 1
fi
WORK=$(mktemp -d)
zcat "$ANNOT22" | awk 'NR>1{print $3}' | sort > "$WORK/annot22.ids"
cut -f2 "$PLINK/1000G.EUR.QC.22.bim" | sort > "$WORK/bim22.ids"
n_annot=$(wc -l < "$WORK/annot22.ids")
n_bim=$(wc -l < "$WORK/bim22.ids")
n_shared=$(comm -12 "$WORK/annot22.ids" "$WORK/bim22.ids" | wc -l)
n_annot_only=$(comm -23 "$WORK/annot22.ids" "$WORK/bim22.ids" | wc -l)
n_bim_only=$(comm -13 "$WORK/annot22.ids" "$WORK/bim22.ids" | wc -l)
frac=$(awk -v a="$n_shared" -v b="$n_bim" 'BEGIN{if(b>0) printf "%.6f", a/b; else print "NA"}')
echo "chr22 baselineLD annot SNPs : $n_annot"
echo "chr22 PLINK panel SNPs      : $n_bim"
echo "chr22 shared identifiers    : $n_shared  (fraction of PLINK covered $frac)"
echo "chr22 annot-only / bim-only : $n_annot_only / $n_bim_only"
rm -f "$WORK/annot22.ids" "$WORK/bim22.ids"; rmdir "$WORK"
if [[ "$n_annot_only" -ne 0 || "$n_bim_only" -ne 0 ]]; then
  echo "chr22 SNP sets are NOT identical: $n_annot_only only in baselineLD, $n_bim_only only in PLINK" >&2
  exit 1
fi
echo "chr22 SNP identifier sets are identical"

cat > "$DEST/DOWNLOAD_COMPLETE" <<EOM
completed $(date -Is) by job ${SLURM_JOB_ID:-unknown} on $(hostname)
1000G_Phase3_baselineLD_v2.2_ldscores.tgz b261e0caf06a003e7522938e01b3d349
1000G_Phase3_frq.tgz                      ac29686ffd5b6378789857a522ebca77
1000G_Phase3_weights_hm3_no_MHC.tgz       a98ac0f089ee285177544a3e6e721ca3
EOM

echo "DOWNLOAD DONE"
