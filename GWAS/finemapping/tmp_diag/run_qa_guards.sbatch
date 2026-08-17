#!/bin/bash
#SBATCH --job-name=colocqa
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/coloc_perf/qaguards_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/coloc_perf/qaguards_%j.err
#
# QA guards against the PROMOTED canonical (2026-08-17). All are READ-ONLY audits
# that exit non-zero on regression. A failure here is the rollback signal;
# rollback is one cp -a from results/susie_coloc_preRerun_backup_2026-08-17.
set +e
set +u
export PATH="${HOME}/.local/bin:${PATH}"
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
QA=$MASLD_PROJECT_ROOT/GWAS/finemapping/results/qa_campaign
cd $MASLD_PROJECT_ROOT/GWAS/finemapping

declare -A RC
run() { # $1 = label, $2.. = command
  local lab="$1"; shift
  echo ""; echo "########## $lab ##########"
  "$@" 2>&1 | tail -40
  RC[$lab]=${PIPESTATUS[0]}
  echo "---- $lab exit=${RC[$lab]} ----"
}
run C4_lambda_s        Rscript $QA/audit/C4_guard.R
run C6_neff            Rscript $QA/audit/C6_neff_audit.R
run C8_controls        Rscript $QA/audit/C8_control_recovery.R
run C9_coverage        Rscript $QA/audit/C9_coverage_audit.R
run C10_reconcile      Rscript $QA/audit/C10_atlas_master_reconcile.R
run provenance         Rscript $QA/generate_provenance_manifest.R
run leakage            bash    $QA/leakage_check.sh

echo ""; echo "================= QA GUARD SUMMARY ================="
fail=0
for k in C4_lambda_s C6_neff C8_controls C9_coverage C10_reconcile provenance leakage; do
  s=$([ "${RC[$k]}" = "0" ] && echo PASS || echo "*** FAIL (exit ${RC[$k]})")
  printf "  %-16s %s\n" "$k" "$s"
  [ "${RC[$k]}" != "0" ] && fail=$((fail+1))
done
echo "  ---------------------------------"
echo "  $((7-fail))/7 passed"
[ "$fail" -gt 0 ] && echo "  *** ROLLBACK SIGNAL: cp -a results/susie_coloc_preRerun_backup_2026-08-17 -> results/susie_coloc"
