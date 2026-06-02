#!/bin/bash -l
# run_recovery_coloc_targeted.sh
# Recover 9 missing loci from the COLOC-targeted fine-mapping run.
#
# Failure modes diagnosed:
#   (a) MHC LD gap in PolyFun panel (6 loci on chr6:28-33 Mb)
#       UKBB_ALT  6.31749142  (array_task=35)
#       UKBB_AST  6.28897980  (array_task=48)
#       UKBB_AST  6.29842451  (array_task=49)
#       UKBB_AST  6.31745464  (array_task=50)
#       UKBB_GGT  6.28898287  (array_task=47)
#       UKBB_GGT  6.31749142  (array_task=48)
#       --> Re-run with LD_PANEL=1kg (1KG EUR fully covers MHC).
#
#   (b) Truncated PolyFun block on chr4 (1 locus)
#       UKBB_ALT  4.52994195  (array_task=26)
#       --> PolyFun chr4 48123600.53877433 only spans 49.0-49.66 Mb on disk;
#           1KG EUR's same-named block properly spans 48.12-53.88 Mb.
#
#   (c) Zero SE in 2023_NAFLD_UKBB_EUR sumstats (2 loci)
#       2023_36280732_NAFLD_UKBB_EUR  2.27598097    (array_task=1)
#       2023_36280732_NAFLD_UKBB_EUR  10.101983413  (array_task=2)
#       --> 16.7% of variants in these loci have se=0; SuSiE rejects.
#           Filter se>0 from existing per-locus ss/ files, then re-run with default LD_PANEL.

set -euo pipefail

FM_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
cd "${FM_DIR}"

mkdir -p logs

# ---------------------------------------------------------------------------
# Step 1: Filter se>0 from 2023_NAFLD_UKBB_EUR per-locus ss/ files (mode c)
# ---------------------------------------------------------------------------
SS_DIR="output/2023_36280732_NAFLD_UKBB_EUR_coloc_targeted/EUR_0.5Mb/ss"
echo "=== Step 1: filtering se>0 from NAFLD_UKBB_EUR ss/ files ==="
for locus in 2.27598097 10.101983413; do
  f="${SS_DIR}/2023_36280732_NAFLD_UKBB_EUR_coloc_targeted_0.5Mb_${locus}.txt"
  if [ ! -f "${f}" ]; then
    echo "  ERROR: ${f} not found, skipping"
    continue
  fi
  # back up once (don't clobber existing backup)
  if [ ! -f "${f}.preSEfilter.bak" ]; then
    cp "${f}" "${f}.preSEfilter.bak"
  fi
  n_before=$(( $(wc -l < "${f}.preSEfilter.bak") - 1 ))
  awk -F'\t' '
    NR == 1 {
      for (i = 1; i <= NF; i++) if ($i == "se") se_col = i
      print
      next
    }
    $se_col + 0 > 0
  ' "${f}.preSEfilter.bak" > "${f}"
  n_after=$(( $(wc -l < "${f}") - 1 ))
  echo "  ${locus}: ${n_after}/${n_before} variants kept ($((n_before - n_after)) dropped for se<=0)"
done

# ---------------------------------------------------------------------------
# Step 2: submit 9 recovery sbatch tasks
# ---------------------------------------------------------------------------
# Common args for run script: <study_name> <ld_pop> <lead_file> <N_tot> <N_cases> <window_mb> <ancestry>

submit_one() {
  local study=$1
  local task_idx=$2
  local ld_panel=$3
  local n_tot=$4
  local n_cases=$5

  local tagged="${study}_coloc_targeted"
  local lead_file="data/lead_snps_coloc_targeted/${study}_leadSNPs.tsv"

  local job_id
  job_id=$(sbatch --parsable \
    --export=ALL,LD_PANEL=${ld_panel} \
    --array=${task_idx} \
    --job-name="recov_${tagged}" \
    --output="logs/recov_${tagged}_%A_%a.out" \
    --error="logs/recov_${tagged}_%A_%a.err" \
    src/03_run_fm_per_locus.sh \
    "${tagged}" "EUR" "${lead_file}" "${n_tot}" "${n_cases}" "0.5" "EUR")
  echo "  ${study} task=${task_idx} LD=${ld_panel} -> job ${job_id}"
}

echo ""
echo "=== Step 2a: 7 LD-recovery tasks (LD_PANEL=1kg) ==="
# UKBB_ALT  N_tot=343850, N_cases=0 (quantitative)
submit_one "UKBB_ALT" 26 "1kg" 343850 0     # chr4:52994195 truncated PolyFun block
submit_one "UKBB_ALT" 35 "1kg" 343850 0     # chr6:31749142 MHC
# UKBB_AST  N_tot=343850, N_cases=0
submit_one "UKBB_AST" 48 "1kg" 343850 0     # chr6:28897980 MHC
submit_one "UKBB_AST" 49 "1kg" 343850 0     # chr6:29842451 MHC
submit_one "UKBB_AST" 50 "1kg" 343850 0     # chr6:31745464 MHC
# UKBB_GGT  N_tot=343850, N_cases=0
submit_one "UKBB_GGT" 47 "1kg" 343850 0     # chr6:28898287 MHC
submit_one "UKBB_GGT" 48 "1kg" 343850 0     # chr6:31749142 MHC

echo ""
echo "=== Step 2b: 2 SE-fix tasks (default LD_PANEL=polyfun, ss files already filtered) ==="
# 2023_36280732_NAFLD_UKBB_EUR  N_tot=400000, N_cases=5765 (binary)
submit_one "2023_36280732_NAFLD_UKBB_EUR" 1 "polyfun" 400000 5765
submit_one "2023_36280732_NAFLD_UKBB_EUR" 2 "polyfun" 400000 5765

echo ""
echo "=== Done. 9 recovery tasks submitted. ==="
echo "Monitor: squeue -u \$(whoami) -o '%.10i %.30j %.2t %.10M' | grep recov_"
