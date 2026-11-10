#!/bin/bash
#SBATCH --job-name=SteatoSITE
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=io,cpu
#SBATCH --prefer=v6
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/%x_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/%x_%j.err
# SS job S2: count assembly (forked r8 assembler, primary arm), identity (s2_identity.py) and blinded QC (s2_qc.R).
# The stage file is not opened by any step. Ends by writing CHECKED_BEFORE_OUTCOME.txt.
set -euo pipefail
SS=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
R=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
PY=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python
# Test mode: S2_TEST=<tag> S2_LIBS=<list> writes everything under verify/s2_test_<tag>/ and no CHECKED file.
libs=${S2_LIBS:-$SS/scripts/s1_samples.txt}
O=$SS; [ -n "${S2_TEST:-}" ] && O=$SS/verify/s2_test_$S2_TEST
for o in out_s2_assembled out_s2_identity out_s2_qc; do [ -e "$O/$o" ] && { echo "refusing to overwrite $O/$o"; exit 1; }; done
mkdir -p "$O"
module load plink/2.0a5.13 >/dev/null 2>&1 && export PLINK2=plink || { module load PLINK/2.00a3.7-gfbf-2023a >/dev/null 2>&1; export PLINK2=plink2; }
sha256sum "$SS"/scripts/s2_*.R "$SS"/scripts/s2_*.py "$SS/scripts/run_s2.sbatch" > "$SS/logs/script_sha256_s2.txt"

nc=0; nf=0; np=0
mkdir -p "$O/out_s2_tmp_$SLURM_JOB_ID"; counted=$O/out_s2_tmp_$SLURM_JOB_ID/counted_libraries.txt; : > "$counted"
while read -r s; do
  if [ -s "$SS/out_s1b_counts/samples/$s/count_completion.json" ] && [ -s "$SS/out_s1c_genotypes/vcf/$s.vcf.gz" ]; then
    echo "$s" >> "$counted"; nc=$((nc+1))
  elif [ -e "$SS/strand_fail_$s" ]; then nf=$((nf+1)); echo "strand fail (excluded before S2): $s"
  elif [ ! -d "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/Human/Patient_Cohorts/results/SteatoSITE/alignments/star/$s" ] && compgen -G "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/Human/Patient_Cohorts/results/SteatoSITE/alignments/star/${s}_ambiguous_*" > /dev/null; then np=$((np+1)); echo "$s" >> "$O/probe_stopped_libraries.txt"; echo "stopped by mate-order probe (excluded before S2): $s"
  else echo "missing: $s" >&2; fi
done < "$libs"
n=$(wc -l < "$libs")
echo "libraries $n counted+genotyped $nc strand_fail $nf probe_stopped $np"
[ $((nc+nf+np)) -eq "$n" ] || { echo "S1 incomplete" >&2; exit 1; }
libs_all=$libs; libs=$counted

roster=$O/out_s2_tmp_$SLURM_JOB_ID/roster.tsv
{ printf 'sample_id\tanalysis_unit_id\tunit_evidence\n'; awk '{print $1"\t"$1"\tsteatosite_library_before_identity_screen"}' "$libs"; } > "$roster"
"$R" --vanilla "$SS/scripts/s2_assemble_verified_counts.R" \
  "$ROOT/figures/candidates/resource-corrections-20261005-v1/annotation_resolution_r8/resolution_summary.json" \
  "$roster" "$SS/out_s1b_counts/samples" "$O/out_s2_assembled" primary
cp "$roster" "$O/out_s2_assembled/roster.tsv"

for s in $(cat "$libs"); do
  v=$SS/out_s1c_genotypes/vcf/$s.vcf.gz
  [ -s "$v.tbi" ] || /gpfs/commons/home/jameslee/micromamba/envs/allelic_qtl/bin/bcftools index -t "$v"
done
"$PY" -I "$SS/scripts/s2_identity.py" --ss-vcf-dir "$SS/out_s1c_genotypes/vcf" \
  --targets-annotated "$SS/out_s0d_targets/s1b_targets_annotated.tsv" \
  --a1-vcf-dir /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/results/rna_genotypes_restricted/a1-full-20260924T025858Z/vcf \
  --libraries "$libs" --out "$O/out_s2_identity" --plink2 "$PLINK2" --threads 8
chmod 750 "$O/out_s2_identity/vcf"
"$PY" -m pip freeze > "$O/out_s2_identity/pip_freeze.txt"

"$R" --vanilla "$SS/scripts/s2_qc.R" "$O/out_s2_assembled" "$O/out_s2_identity" "$O/out_s2_qc"
rm -rf "$O/out_s2_tmp_$SLURM_JOB_ID"

[ -n "${S2_TEST:-}" ] && { echo "S2 test done: $O"; exit 0; }
{ echo "SS S2 blinded QC completed $(date -u +%Y-%m-%dT%H:%M:%SZ), job $SLURM_JOB_ID. The stage file was not opened."
  echo "Scripts: $(cat "$SS/logs/script_sha256_s2.txt" | awk '{print $2": "$1}' | tr '\n' ';')"
  echo "units.tsv sha256: $(sha256sum "$O/out_s2_qc/units.tsv" | cut -c1-64)"
  echo "dge_included.rds sha256: $(sha256sum "$O/out_s2_qc/dge_included.rds" | cut -c1-64)"
  echo "Libraries listed: $n; excluded before S2 for read 1 antisense fraction <= 0.80: $nf"
  cat "$O/out_s2_qc/S2_summary.json"; } > "$SS/CHECKED_BEFORE_OUTCOME.txt"
echo "S2 done"
