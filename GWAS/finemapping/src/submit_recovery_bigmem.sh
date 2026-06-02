#!/bin/bash
# submit_recovery_bigmem.sh — submits 4 bigmem interactive-QOS jobs
# for SuSiE-COLOC recovery. Auto-detects bigmem free mem and picks tier.

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping

# --- Detect bigmem availability ---
free_mem_gb=$(sinfo -p bigmem -h -N -o "%e" | awk '{s+=$1} END {printf "%d", s/1024}')
echo "Bigmem free memory across nodes: ${free_mem_gb} GB"

# --- Pick memory tier (6 TB cap) ---
if [ "$free_mem_gb" -ge 5500 ]; then
  MEM_PER_JOB=1500; WORKERS_SAFE=32; WORKERS_RISKY=16; TIER="full"
elif [ "$free_mem_gb" -ge 4200 ]; then
  MEM_PER_JOB=1100; WORKERS_SAFE=24; WORKERS_RISKY=14; TIER="moderate"
elif [ "$free_mem_gb" -ge 3000 ]; then
  MEM_PER_JOB=800;  WORKERS_SAFE=20; WORKERS_RISKY=12; TIER="conservative"
else
  MEM_PER_JOB=600;  WORKERS_SAFE=16; WORKERS_RISKY=10; TIER="minimal"
fi
total_mb=$((MEM_PER_JOB * 4))
echo "Selected tier: ${TIER} (${MEM_PER_JOB}G/job × 4 = ${total_mb}G total)"
echo "Workers: safe=${WORKERS_SAFE} risky=${WORKERS_RISKY}"
echo ""

# --- Submit 4 jobs ---
declare -A TASKS_WORKERS=(
  ["data/recovery/tasks_a_safe.txt"]="${WORKERS_SAFE}"
  ["data/recovery/tasks_b_safe.txt"]="${WORKERS_SAFE}"
  ["data/recovery/tasks_a_risky.txt"]="${WORKERS_RISKY}"
  ["data/recovery/tasks_b_risky.txt"]="${WORKERS_RISKY}"
)

mkdir -p logs
echo "Submitted jobs:"
for task_file in data/recovery/tasks_a_safe.txt data/recovery/tasks_b_safe.txt \
                  data/recovery/tasks_a_risky.txt data/recovery/tasks_b_risky.txt; do
  workers=${TASKS_WORKERS[$task_file]}
  n_tasks=$(wc -l < "$task_file")
  jid=$(TASK_FILE="$task_file" N_WORKERS="$workers" \
    sbatch --parsable \
      --mem=${MEM_PER_JOB}G \
      --job-name=recov_$(basename "$task_file" .txt) \
      src/06_susie_coloc_bigmem_recovery.sh)
  echo "  ${jid}  ${task_file}  (${n_tasks} tasks, ${workers} workers, ${MEM_PER_JOB}G)"
done

echo ""
echo "Monitor: squeue -u \$(whoami) -n recov_tasks_a_safe,recov_tasks_b_safe,recov_tasks_a_risky,recov_tasks_b_risky"
echo "Or:     squeue -u \$(whoami) | grep recov"
