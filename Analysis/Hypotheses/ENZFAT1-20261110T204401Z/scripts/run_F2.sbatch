#!/bin/bash
#SBATCH --job-name=ENZFAT1_F2
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --time=90:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/ENZFAT1-20261110T204401Z/logs/F2_%j.out
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
RUN=$ROOT/Analysis/Hypotheses/ENZFAT1-20261110T204401Z
OUT=$RUN/out_F2
S3=$RUN/out_S3_r2
PY=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python
LD=/gpfs/commons/home/jameslee/micromamba/envs/ldsc/bin
REF=$ROOT/GWAS/ldsc_ref/eur_w_ld_chr
MA=$ROOT/GWAS/ldsc_ref/hm3_merge_alleles.txt
TRAITS="bmi tg fi dpw t2d height"
if [[ -e $OUT ]]; then echo "$OUT exists; not overwriting"; exit 1; fi
mkdir -p $OUT/inputs $OUT/checks $OUT/munged $OUT/ldsc
(cd $RUN && head -1 PRESPEC_ENZFAT1_ADDENDUM2_F2.sha256 | sha256sum -c -)
sha256sum $RUN/PRESPEC_ENZFAT1_ADDENDUM2_F2.md $RUN/scripts/40_build_f2.py $RUN/scripts/41_f2_checks.py $RUN/scripts/42_f2_summary.py $RUN/scripts/run_F2.sbatch $S3/ldsc/component.sumstats.gz > $OUT/script_and_input_sha256.txt
$PY -I -m pip freeze > $OUT/pip_freeze_spatial.txt
echo "seed 20261110 (no random step)" > $OUT/seed.txt
# S3 gate (Addendum 2 section 2): component h2 z >= 4 and chr22 userGWAS z correlation > 0.99
$PY -I - <<PYEOF
import re, sys, pandas as pd
t = open("$S3/ldsc/component_h2.log").read()
h, s = map(float, re.search(r"Total Observed scale h2: ([-\d.e]+) \(([-\d.e]+)\)", t).groups())
c = pd.read_csv("$S3/usergwas_chr22_check.tsv", sep="\t").iloc[0]
ok = h / s >= 4 and bool(c["pass_cor_gt_0.99"])
print(f"S3 gate: component h2 z {h/s:.2f}; chr22 cor {c['cor_z']:.4f}; gate {'pass' if ok else 'FAIL'}")
sys.exit(0 if ok else 1)
PYEOF
# md5 checks
(cd $RUN/downloads_untrusted/tg_GCST90239664 && grep buildGRCh37.tsv$ md5sum.txt | md5sum -c -)
(cd $RUN/downloads_untrusted/height_GCST90245992 && grep buildGRCh37.tsv$ md5sum.txt | md5sum -c -)
(cd $RUN/downloads_untrusted/fi_GCST90002238 && grep buildGRCh37.tsv.gz$ md5sum.txt | md5sum -c -)
echo "b19a25d153470f0de906a2ac888a936e  $ROOT/Analysis/Hypotheses/DIR6-20261111T024808Z/downloads_untrusted/GCST90475667/GCST90475667.h.tsv.gz" | md5sum -c -
md5sum $RUN/downloads_untrusted/dpw_GCST007461/DrinksPerWeek.txt.gz > $OUT/dpw_local_md5.txt
echo $TRAITS | tr ' ' '\n' | xargs -P 4 -I{} sh -c "$PY -I $RUN/scripts/40_build_f2.py {} $ROOT $RUN $OUT > $OUT/inputs/{}.log 2>&1"
cat $OUT/inputs/*.log
$PY -I $RUN/scripts/41_f2_checks.py $OUT
for t in $TRAITS; do
  $LD/python $LD/munge_sumstats.py --sumstats $OUT/inputs/$t.tsv.gz --snp SNP --a1 A1 --a2 A2 \
    --signed-sumstats BETA,0 --p P --N-col N --merge-alleles $MA --chunksize 500000 --out $OUT/munged/$t \
    > $OUT/munged/$t.stdout 2>&1
  $LD/python $LD/ldsc.py --h2 $OUT/munged/$t.sumstats.gz --ref-ld-chr $REF/ --w-ld-chr $REF/ --out $OUT/ldsc/h2_$t > /dev/null 2>&1
done
L=""; for t in $TRAITS; do L="$L,$OUT/munged/$t.sumstats.gz"; done
$LD/python $LD/ldsc.py --rg $S3/ldsc/component.sumstats.gz$L --ref-ld-chr $REF/ --w-ld-chr $REF/ --out $OUT/ldsc/rg_component > /dev/null 2>&1
$LD/python $LD/ldsc.py --rg $RUN/out_S0/munged/alt_ukb.sumstats.gz$L --ref-ld-chr $REF/ --w-ld-chr $REF/ --out $OUT/ldsc/rg_alt_ukb > /dev/null 2>&1
$LD/python $LD/ldsc.py --rg $RUN/out_S0/munged/pdff.sumstats.gz$L --ref-ld-chr $REF/ --w-ld-chr $REF/ --out $OUT/ldsc/rg_pdff > /dev/null 2>&1
for t in $TRAITS; do echo "$t $(grep 'Total Observed scale h2' $OUT/ldsc/h2_$t.log)"; done
$PY -I $RUN/scripts/42_f2_summary.py $OUT
echo "F2 DONE"
