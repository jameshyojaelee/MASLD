#!/bin/bash
# Runs the post-completion chain the moment the COLOC rerun is genuinely done.
#
# COMPLETION IS OUTPUT-BASED, NOT STATE-BASED. 42 tasks are permanently FAILED
# (40 died after writing their results during a 35-min script corruption on
# Aug 11; 2 were repaired and their outputs landed), so "no jobs left" and "all
# work done" are different questions. This waits for BOTH: the queue to drain
# AND 1100 output files to exist. If the queue empties while outputs are short,
# it says so and stops rather than proceeding on incomplete data -- pre-
# registration rule 5.
#
# Chain, in order, each gated on the previous:
#   1. structural verification of all 1100 outputs   (read-only)
#   2. pre-registered comparison P1-P5 + P7          (read-only)
#   3. P6 via the REDIRECTED copy of 07              (writes only into _rerun)
#
# Step 3 uses src/perf/07_combine_rerun.R, which carries a refuse-to-run guard:
# the original hardcodes COLOC_DIR=results/susie_coloc and would overwrite
# gene_level_coloc.csv, the firewall-guarded file behind the atlas and Fig2.
FM=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
LOG="$FM/tmp_diag/completion_chain.log"
log() { echo "$(date -u +%FT%TZ) $*" | tee -a "$LOG"; }

log "waiting for completion (queue drained AND 1100 outputs)"
while :; do
  q=$(squeue -u "$USER" -h -r -n susie_coloc -t RUNNING,PENDING 2>/dev/null | wc -l)
  n=$(find "$FM/results/susie_coloc_rerun" -name 'susie_coloc_chr*.csv' 2>/dev/null | wc -l)
  if [ "$q" -eq 0 ] && [ "$n" -ge 1100 ]; then log "COMPLETE: queue empty, $n/1100 outputs"; break; fi
  if [ "$q" -eq 0 ] && [ "$n" -lt 1100 ]; then
    log "HALT: queue empty but only $n/1100 outputs. The watcher should repair;"
    log "      holding rather than analysing incomplete data (prereg rule 5)."
    sleep 600; continue
  fi
  sleep 600
done

cd "$FM"

log "STEP 1/3  structural verification"
J1=$(sbatch --parsable src/perf/verify_rerun_outputs.sbatch)
while squeue -j "$J1" -h -t RUNNING,PENDING 2>/dev/null | grep -q .; do sleep 60; done
log "  job $J1 done -> logs/coloc_perf/verify_outputs_${J1}.out"
if grep -q "PROBLEMS FOUND" "logs/coloc_perf/verify_outputs_${J1}.out" 2>/dev/null; then
  log "  *** STRUCTURAL PROBLEMS -- stopping the chain; do not analyse these outputs"
  exit 1
fi

log "STEP 2/3  pre-registered comparison (P1-P5, P7)"
J2=$(sbatch --parsable src/perf/prereg_comparison.sbatch)
while squeue -j "$J2" -h -t RUNNING,PENDING 2>/dev/null | grep -q .; do sleep 60; done
log "  job $J2 done -> logs/coloc_perf/prereg_${J2}.out"

log "STEP 3/3  P6: gene-level tier table from the rerun (redirected copy of 07)"
cat > /tmp/p6_$$.sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=coloccmb
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=$FM/logs/coloc_perf/p6_%j.out
#SBATCH --error=$FM/logs/coloc_perf/p6_%j.err
set -eo pipefail
set +u
export PATH="\${HOME}/.local/bin:\${PATH}"
eval "\$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -u
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd $FM/src
Rscript $FM/src/perf/07_combine_rerun.R
EOF
J3=$(sbatch --parsable /tmp/p6_$$.sbatch); rm -f /tmp/p6_$$.sbatch
while squeue -j "$J3" -h -t RUNNING,PENDING 2>/dev/null | grep -q .; do sleep 60; done
log "  job $J3 done -> logs/coloc_perf/p6_${J3}.out"

log "CHAIN COMPLETE"
log "  verification : logs/coloc_perf/verify_outputs_${J1}.out"
log "  comparison   : logs/coloc_perf/prereg_${J2}.out"
log "  P6 table     : logs/coloc_perf/p6_${J3}.out"
log "  canonical gene_level_coloc.csv untouched: $(ls -l $FM/results/susie_coloc/gene_level_coloc.csv 2>/dev/null | awk '{print $6,$7,$8}')"
