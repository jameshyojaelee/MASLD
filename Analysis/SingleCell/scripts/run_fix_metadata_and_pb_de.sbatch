#!/bin/bash
#SBATCH --job-name=fix_meta_pb
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/logs/fix_meta_pb_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/logs/fix_meta_pb_%j.err

set -euo pipefail
BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "$BASE"
mkdir -p Analysis/SingleCell/logs

echo "=== Step 1a: Add condition_harmonized + preparation_method to atlas ==="
echo "Start: $(date)"

module purge 2>/dev/null || true
eval "$(micromamba shell hook -s bash)"
micromamba activate rapids_singlecell
export PYTHONNOUSERSITE=1

python Analysis/SingleCell/integration/scripts/add_sample_metadata.py

echo ""
echo "=== Verify metadata columns ==="
python -c "
import h5py
f = h5py.File('Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad', 'r')
obs_keys = list(f['obs'].keys())
print('obs columns:', obs_keys)
assert 'condition_harmonized' in obs_keys, 'condition_harmonized MISSING'
assert 'preparation_method' in obs_keys, 'preparation_method MISSING'
print('OK: Both columns present')
f.close()
"

echo ""
echo "=== Step 1b: Re-run pseudobulk DE with fixed metadata ==="
echo "Start: $(date)"

micromamba activate rnaseq
Rscript Analysis/SingleCell/scripts/pseudobulk_de.R

echo ""
echo "=== Step 1c: Evaluate new correlation ==="
micromamba activate rapids_singlecell
export PYTHONNOUSERSITE=1

python3 -c "
import pandas as pd, numpy as np, os

# Gene map for symbol->ensembl
atlas = pd.read_csv('RNA-seq/results/multi_evidence/multi_evidence_atlas.csv', usecols=['ensembl_id','human_symbol'])
atlas['ensembl_clean'] = atlas['ensembl_id'].str.replace(r'\..*', '', regex=True)
sym_map = atlas.dropna(subset=['human_symbol']).drop_duplicates(subset='human_symbol')[['human_symbol','ensembl_clean']]

dream = pd.read_csv('RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/dream_results.csv')
dream = dream.rename(columns={'logFC':'dream_logFC', 'padj':'dream_padj'})
dream['ensembl_clean'] = dream['gene'].str.replace(r'\..*', '', regex=True)

pb_dir = 'Analysis/SingleCell/results_gpu_v2/pseudobulk_de'
print(f'{'Cell type':<30} {'r':>6} {'n_genes':>8} {'both_sig':>9} {'pb_sig':>7}')
print('-' * 65)
for f in sorted(os.listdir(pb_dir)):
    if not f.endswith('_de.csv'): continue
    ct = f.replace('_de.csv','')
    pb = pd.read_csv(os.path.join(pb_dir, f))
    is_ensg = pb['gene'].str.startswith('ENSG')
    pb.loc[is_ensg, 'ensembl_clean'] = pb.loc[is_ensg, 'gene'].str.replace(r'\..*', '', regex=True)
    pb = pb.merge(sym_map, left_on='gene', right_on='human_symbol', how='left', suffixes=('','_mapped'))
    pb.loc[~is_ensg & pb['ensembl_clean_mapped'].notna(), 'ensembl_clean'] = pb.loc[~is_ensg & pb['ensembl_clean_mapped'].notna(), 'ensembl_clean_mapped']
    pb = pb.dropna(subset=['ensembl_clean'])
    m = pb.merge(dream[['ensembl_clean','dream_logFC','dream_padj']], on='ensembl_clean')
    valid = m.dropna(subset=['logFC','dream_logFC'])
    if len(valid) < 100:
        print(f'{ct:<30} SKIP (n={len(valid)})')
        continue
    r = np.corrcoef(valid['logFC'], valid['dream_logFC'])[0,1]
    both = ((valid['padj'] < 0.05) & (valid['dream_padj'] < 0.1)).sum()
    pb_sig = (valid['padj'] < 0.05).sum()
    print(f'{ct:<30} {r:>6.3f} {len(valid):>8} {both:>9} {pb_sig:>7}')
"

echo ""
echo "=== Done: $(date) ==="
