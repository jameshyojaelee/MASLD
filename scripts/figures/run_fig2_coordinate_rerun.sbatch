#!/usr/bin/env bash
# Corrected-coordinate candidates only; never overwrite a named release.
#SBATCH --job-name=fig2_coords
#SBATCH --account=nslab
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/logs/fig2_coordinates_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/logs/fig2_coordinates_%j.err

set -euo pipefail
umask 027
project_root=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
python_bin=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
run_stamp=$(date -u +%Y%m%dT%H%M%SZ)
run_root="${project_root}/figures/candidates/fig2-coordinate-rerun-${run_stamp}_${SLURM_JOB_ID}"
snapshot_root="${run_root}/snapshot"
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export FIG2_PROJECT_ROOT="${project_root}"
test ! -e "${run_root}"
mkdir "${run_root}"
mkdir -p "${snapshot_root}/scripts/figures/tests"
printf '%s\n' "${run_root}" > "${project_root}/logs/fig2_coordinates_${SLURM_JOB_ID}.run_root.txt"
printf '[fig2] job=%s host=%s run=%s start=%s\n' "${SLURM_JOB_ID}" "$(hostname)" "${run_root}" "$(date -u)"

for script_name in fig2_locus_coordinates.py test_fig2_locus_coordinates.py \
    fig2_multiplicity_correction.py fig2_coding_mass_robustness.py \
    fig2_trait_directness_stats.py fig2_ascertainment_control.py; do
    cp -p "${project_root}/scripts/figures/${script_name}" "${snapshot_root}/scripts/figures/"
done
cp -p "${project_root}/scripts/figures/tests/test_fig2_cluster_inference.py" \
    "${snapshot_root}/scripts/figures/tests/"
cp -p "${project_root}/scripts/figures/run_fig2_coordinate_rerun.sbatch" "${run_root}/launcher.sbatch"
"${python_bin}" -B -c 'import numpy, pandas, statsmodels; import sys; print(sys.version); print("numpy", numpy.__version__); print("pandas", pandas.__version__); print("statsmodels", statsmodels.__version__)' > "${run_root}/environment.txt"
"${python_bin}" -m pip freeze >> "${run_root}/environment.txt"

# Small source tables are copied for the record. Large summary statistics are
# read in place and hashed before and after execution on this compute node.
"${python_bin}" -B - "${project_root}" "${run_root}" <<'PY'
import csv
import hashlib
import pathlib
import shutil
import sys

root, run = map(pathlib.Path, sys.argv[1:])
base = pathlib.Path('GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321')
registry_rel = pathlib.Path('GWAS/finemapping/config/gwas_registry.tsv')
small = [base / name for name in ('credible_set_pip_architecture.tsv',
         'credible_set_members_annotated.tsv', 'credible_set_regulatory_context.tsv')]
small.append(registry_rel)
small.extend(pathlib.Path('GWAS/finemapping/runs/uniform35_v3_2026-08-03/config/ld_blocks') / f'{anc}.tsv'
             for anc in ('AFR', 'AMR', 'EAS', 'EUR', 'SAS'))
for rel in small:
    target = run / 'snapshot' / rel
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(root / rel, target)
with (root / base / 'credible_set_pip_architecture.tsv').open() as handle:
    studies = {row['study'] for row in csv.DictReader(handle, delimiter='\t')
               if row['trait_scope'] in ('tier1_direct_masld_pdff', 'tier2_liver_enzyme')}
with (root / registry_rel).open() as handle:
    registry = list(csv.DictReader(handle, delimiter='\t'))
sumstats = {root / 'GWAS/finemapping/data/sumstats' / pathlib.Path(row['sumstats_path']).name
            for row in registry if row['study_name'] in studies}
paths = {root / rel for rel in small} | sumstats
paths |= {path for path in (run / 'snapshot').rglob('*') if path.is_file()}
paths.add(run / 'launcher.sbatch')
rows = []
for path in sorted(paths):
    if not path.exists():
        rows.append({'path': str(path), 'size_bytes': '', 'sha256': '', 'state': 'missing'})
        continue
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b''):
            digest.update(chunk)
    rows.append({'path': str(path), 'size_bytes': path.stat().st_size,
                 'sha256': digest.hexdigest(), 'state': 'present'})
with (run / 'input_checksums_before.tsv').open('w') as handle:
    writer = csv.DictWriter(handle, fieldnames=('path', 'size_bytes', 'sha256', 'state'), delimiter='\t')
    writer.writeheader()
    writer.writerows(rows)
PY

printf '[fig2] focused coordinate and cluster tests start=%s\n' "$(date -u)"
"${python_bin}" -B "${snapshot_root}/scripts/figures/test_fig2_locus_coordinates.py" \
    > "${run_root}/coordinate_tests.log" 2>&1
"${python_bin}" -B "${snapshot_root}/scripts/figures/tests/test_fig2_cluster_inference.py" \
    > "${run_root}/cluster_tests.log" 2>&1

for analysis_name in multiplicity robustness trait_directness ascertainment; do
    case "${analysis_name}" in
        multiplicity) producer=fig2_multiplicity_correction.py ;;
        robustness) producer=fig2_coding_mass_robustness.py ;;
        trait_directness) producer=fig2_trait_directness_stats.py ;;
        ascertainment) producer=fig2_ascertainment_control.py ;;
    esac
    printf '[fig2] analysis=%s start=%s\n' "${analysis_name}" "$(date -u)"
    "${python_bin}" -B "${snapshot_root}/scripts/figures/${producer}" \
        --project-root "${project_root}" --output-dir "${run_root}/${analysis_name}" \
        > "${run_root}/${analysis_name}.log" 2>&1
    printf '[fig2] analysis=%s finished=%s\n' "${analysis_name}" "$(date -u)"
done
test -s "${run_root}/multiplicity/corrected_locus_coordinates.tsv"
test -s "${run_root}/robustness/corrected_locus_coordinates.tsv"

"${python_bin}" -B - "${run_root}" <<'PY'
import csv
import hashlib
import pathlib
import sys

run = pathlib.Path(sys.argv[1])
with (run / 'input_checksums_before.tsv').open() as handle:
    rows = list(csv.DictReader(handle, delimiter='\t'))
changed = []
after = []
for row in rows:
    path = pathlib.Path(row['path'])
    if not path.exists():
        current = {'path': str(path), 'size_bytes': '', 'sha256': '', 'state': 'missing'}
    else:
        digest = hashlib.sha256()
        with path.open('rb') as handle:
            for chunk in iter(lambda: handle.read(1 << 20), b''):
                digest.update(chunk)
        current = {'path': str(path), 'size_bytes': str(path.stat().st_size),
                   'sha256': digest.hexdigest(), 'state': 'present'}
    after.append(current)
    if current != row:
        changed.append(str(path))
with (run / 'input_checksums_after.tsv').open('w') as handle:
    writer = csv.DictWriter(handle, fieldnames=('path', 'size_bytes', 'sha256', 'state'), delimiter='\t')
    writer.writeheader()
    writer.writerows(after)
if changed:
    raise RuntimeError('Inputs changed during execution: ' + ', '.join(changed))
(run / 'COMPLETE.txt').write_text('Four corrected-coordinate candidate analyses completed; input checksums unchanged.\n')
PY
printf '[fig2] COMPLETE run=%s finished=%s\n' "${run_root}" "$(date -u)"
