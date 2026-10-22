#!/bin/bash
#SBATCH --job-name=ldsc_pdff_fix
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --time=90:00:00
#SBATCH --output=GWAS/ldsc/rg_pdff_corrected_20261012T012936Z/slurm_%j.out
#SBATCH --error=GWAS/ldsc/rg_pdff_corrected_20261012T012936Z/slurm_%j.err
# Rerun of every LDSC call in GWAS/ldsc/rg/ that reads pdff.sumstats.gz, after the
# 2026-10-10 allele-label correction of the Pazoki 2022 PDFF GWAS (GCST90267352).
# Commands are copied from the "Call:" block of each old log; only --out changes.
# LDSC 1.0.1, micromamba env "ldsc" (python 2.7), 1000G EUR LD scores (eur_w_ld_chr).
set -uo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
OUT=GWAS/ldsc/rg_pdff_corrected_20261012T012936Z
REF=GWAS/ldsc_ref/eur_w_ld_chr; MG=GWAS/ldsc/munged
export MAMBA_ROOT_PREFIX=/gpfs/commons/home/jameslee/micromamba
run(){ /gpfs/commons/home/jameslee/.local/bin/micromamba run -n ldsc "$@"; }

# refuse to run on anything but the corrected munged file
want=0d2e4e00efedd6207ef759068d377b23
got=$(md5sum $MG/pdff.sumstats.gz | cut -d' ' -f1)
[[ "$got" == "$want" ]] || { echo "pdff.sumstats.gz md5 $got != corrected $want; abort"; exit 1; }
md5sum $MG/*.sumstats.gz $MG/pdff.sumstats.gz.as_downloaded_20261010T024403Z > $OUT/input_md5.txt
run pip freeze > $OUT/pip_freeze.txt 2>&1
run python --version > $OUT/python_version.txt 2>&1

s(){ local p=""; for t in "$@"; do p+="${p:+,}$MG/$t.sumstats.gz"; done; echo "$p"; }
rg(){ local out=$1; shift
  run ldsc.py --rg "$(s "$@")" --ref-ld-chr $REF/ --w-ld-chr $REF/ --out $OUT/$out > $OUT/$out.stdout 2>&1 \
    && echo "rg $out ok" || echo "rg $out FAIL"; }

run ldsc.py --h2 $MG/pdff.sumstats.gz --ref-ld-chr $REF/ --w-ld-chr $REF/ --out $OUT/h2_pdff > $OUT/h2_pdff.stdout 2>&1 && echo "h2 ok" || echo "h2 FAIL"

# run_ldsc3.sh set (logs ma_*)
rg ma_ghodsian_nafld ghodsian_nafld mvp_nafld finngen_nafld pdff ghouse_cirr mvp_cirr finngen_cirr ukbb_alt mvp_alt ukbb_ast mvp_ast finngen_obesity
rg ma_mvp_nafld      mvp_nafld finngen_nafld pdff ghouse_cirr mvp_cirr finngen_cirr ukbb_alt mvp_alt ukbb_ast mvp_ast finngen_obesity
rg ma_finngen_nafld  finngen_nafld pdff ghouse_cirr mvp_cirr finngen_cirr ukbb_alt mvp_alt ukbb_ast mvp_ast finngen_obesity
rg ma_pdff           pdff ghouse_cirr mvp_cirr finngen_cirr ukbb_alt mvp_alt ukbb_ast mvp_ast finngen_obesity
# run_ldsc2.sh set (logs wp_*). The July wp_* runs read munged files made without
# --merge-alleles; those were deleted by run_ldsc3.sh, so these reruns read the
# current (--merge-alleles) munged files.
rg wp_ghodsian_nafld ghodsian_nafld mvp_nafld pdff ghouse_cirr mvp_cirr ukbb_alt mvp_alt ukbb_ast mvp_ast
rg wp_mvp_nafld      mvp_nafld pdff ghouse_cirr mvp_cirr ukbb_alt mvp_alt ukbb_ast mvp_ast
rg wp_pdff           pdff ghouse_cirr mvp_cirr ukbb_alt mvp_alt ukbb_ast mvp_ast
echo "DONE"
