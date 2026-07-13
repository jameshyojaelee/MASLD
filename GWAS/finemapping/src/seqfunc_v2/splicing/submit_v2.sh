#!/usr/bin/env bash
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SRC=${ROOT}/GWAS/finemapping/src/seqfunc_v2/splicing
OUT=${ROOT}/GWAS/finemapping/results/seqfunc/disease_splicing/joint_v2
export MASLD_PROJECT_ROOT=${ROOT}
micromamba run -n rnaseq Rscript "${SRC}/00_prepare_manifest.R"
test "$(($(wc -l < "${OUT}/canonical_846_manifest.tsv")-1))" -eq 846

cluster=$(sbatch --parsable "${SRC}/01_joint_cluster.sbatch")
cohort=$(sbatch --parsable "${SRC}/04_cohort_native.sbatch")
primary=$(sbatch --parsable --dependency=afterok:${cluster} "${SRC}/02_primary_ds.sbatch")
loco=$(sbatch --parsable --dependency=afterok:${cluster} "${SRC}/03_loco_ds.sbatch")
perm=$(sbatch --parsable --dependency=afterok:${cluster} "${SRC}/05_permutation_ds.sbatch")
perm_gate=$(sbatch --parsable --dependency=afterok:${primary}:${perm} "${SRC}/06_permutation_gate.sbatch")
annotate=$(sbatch --parsable --dependency=afterok:${primary}:${loco}:${cohort} "${SRC}/07_annotate_replicate.sbatch")
release=$(sbatch --parsable --dependency=afterok:${annotate}:${perm_gate} "${SRC}/08_release_gate.sbatch")
printf 'stage\tjob_id\tdependency\ncluster\t%s\tnone\ncohort\t%s\tnone\nprimary\t%s\tafterok:%s\nloco\t%s\tafterok:%s\npermutation10\t%s\tafterok:%s\npermutation_gate\t%s\tafterok:%s:%s\nannotation\t%s\tafterok:%s:%s:%s\nrelease\t%s\tafterok:%s:%s\n' \
 "$cluster" "$cohort" "$primary" "$cluster" "$loco" "$cluster" "$perm" "$cluster" \
 "$perm_gate" "$primary" "$perm" "$annotate" "$primary" "$loco" "$cohort" "$release" "$annotate" "$perm_gate" \
 | tee "${OUT}/submitted_jobs.tsv"
