#!/bin/bash
#SBATCH --job-name=lncstrand
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=72:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/noncoding_resource/slurm_logs/lncstrand_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/noncoding_resource/slurm_logs/lncstrand_%j.err

set -eo pipefail

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
CANDIDATE_ROOT="${PROJECT_ROOT}/RNA-seq/results/noncoding_resource/candidates/noncoding-bulk-2026-08-11"
RUN_ROOT="${CANDIDATE_ROOT}/strand_audit/GSE162694/run_${SLURM_JOB_ID}"
GTF="${PROJECT_ROOT}/RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BG001-DECISION/external_reference/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
BAM_ROOT="${PROJECT_ROOT}/RNA-seq/Human/Patient_Cohorts/results/GSE162694/alignments/star"

if [[ -e "${RUN_ROOT}" ]]; then
    echo "Refusing to overwrite existing run directory: ${RUN_ROOT}" >&2
    exit 1
fi
mkdir -p "${RUN_ROOT}"

eval "$(micromamba shell hook --shell=bash)"
micromamba activate rnaseq
set -u

version="$(featureCounts -v 2>&1 | awk 'NF {print $NF; exit}')"
if [[ "${version}" != "v2.1.1" ]]; then
    echo "Expected featureCounts v2.1.1, found ${version}" >&2
    exit 1
fi

samples=(
    SRR13199852
    SRR13199853
    SRR13199854
    SRR13199855
    SRR13199856
    SRR13199857
    SRR13199858
    SRR13199859
    SRR13199860
    SRR13199862
)

bams=()
for sample in "${samples[@]}"; do
    bam="${BAM_ROOT}/${sample}/${sample}.Aligned.sortedByCoord.out.bam"
    if [[ ! -s "${bam}" ]]; then
        echo "Missing BAM: ${bam}" >&2
        exit 1
    fi
    bams+=("${bam}")
done

{
    printf 'sample_id\tbam_path\tbam_bytes\n'
    for i in "${!samples[@]}"; do
        printf '%s\t%s\t%s\n' "${samples[$i]}" "${bams[$i]}" "$(stat -c '%s' "${bams[$i]}")"
    done
} > "${RUN_ROOT}/sample_manifest.tsv"

{
    printf 'slurm_job_id\t%s\n' "${SLURM_JOB_ID}"
    printf 'started_utc\t%s\n' "$(date -u +%FT%TZ)"
    printf 'hostname\t%s\n' "$(hostname)"
    printf 'featureCounts_version\t%s\n' "${version}"
    printf 'featureCounts_path\t%s\n' "$(command -v featureCounts)"
    printf 'gtf_path\t%s\n' "${GTF}"
    printf 'gtf_sha256\t%s\n' "$(sha256sum "${GTF}" | awk '{print $1}')"
} > "${RUN_ROOT}/environment.tsv"

for strand in 0 1 2; do
    out="${RUN_ROOT}/gene_counts_s${strand}.txt"
    featureCounts \
        -T "${SLURM_CPUS_PER_TASK}" \
        -s "${strand}" \
        -a "${GTF}" \
        -o "${out}" \
        "${bams[@]}" \
        > "${RUN_ROOT}/featureCounts_s${strand}.stdout.log" \
        2> "${RUN_ROOT}/featureCounts_s${strand}.stderr.log"
done

sha256sum "${RUN_ROOT}"/gene_counts_s*.txt* > "${RUN_ROOT}/outputs.sha256"
printf 'completed_utc\t%s\n' "$(date -u +%FT%TZ)" >> "${RUN_ROOT}/environment.tsv"
touch "${RUN_ROOT}/COMPLETE"
