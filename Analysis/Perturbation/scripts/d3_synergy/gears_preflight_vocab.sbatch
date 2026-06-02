#!/bin/bash
#SBATCH --job-name=GEARS
#SBATCH --partition=cpu
#SBATCH --time=48:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --output=Analysis/Perturbation/scripts/d3_synergy/logs/gears_preflight_vocab_%j.log

# Pre-flight: load Norman processed AnnData (~5.6 GB), extract pert vocab,
# check overlap with D3 tier1/tier2 invocab hits + Saunders panel.
# Gate per team-lead spec: if tier1 overlap ≥3/4 → green light for fine-tune;
# if <2 → halt-and-report deeper mismatch.

set -o pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export PYTHONNOUSERSITE=1
unset PYTHONPATH || true

eval "$(micromamba shell hook --shell bash)"
micromamba activate -p .envs/perturbation_gears

set -u

python - <<'PY'
import anndata as ad
import pandas as pd
import json
from pathlib import Path

NORMAN_H5AD = Path("Analysis/Perturbation/results/finetuned_checkpoints/norman/perturb_processed.h5ad")
print(f"Loading Norman processed AnnData from {NORMAN_H5AD}")
print(f"  size: {NORMAN_H5AD.stat().st_size/1e9:.2f} GB")
adata = ad.read_h5ad(NORMAN_H5AD)
print(f"Norman: {adata.shape}")
print(f"  obs cols: {adata.obs.columns.tolist()}")
print(f"  var cols: {adata.var.columns.tolist()[:10]}")

# Extract pert vocab from condition column
conds = adata.obs['condition'].value_counts()
print(f"\nNorman total unique conditions: {len(conds)}")
print(f"Top 10 conditions:\n{conds.head(10)}")

pert_genes = set()
for c in adata.obs['condition'].unique():
    if c == 'ctrl': continue
    for p in c.split('+'):
        if p != 'ctrl':
            pert_genes.add(p)
print(f"\nUnique Norman perturbed genes (single + combo combined): {len(pert_genes)}")
print(f"Sample 30 sorted: {sorted(pert_genes)[:30]}")

# Tier1 overlap
t1_genes = ['TADA1','VPS41','HECTD1','RNF103']
overlap_t1 = sorted(set(t1_genes) & pert_genes)
print(f"\n=== TIER1 OVERLAP === {len(overlap_t1)}/4 → {overlap_t1}")

# Tier2 overlap
t2 = pd.read_csv('Analysis/Perturbation/data/hits/d3_tier2_pairs_invocab.csv')
t2_genes = sorted(set(t2['gene1']) | set(t2['gene2']))
overlap_t2 = sorted(set(t2_genes) & pert_genes)
t2_predictable = t2[t2['gene1'].isin(pert_genes) & t2['gene2'].isin(pert_genes)]
print(f"=== TIER2 OVERLAP === unique genes: {len(t2_genes)} total, {len(overlap_t2)} in Norman")
print(f"Tier2 PAIRS predictable (BOTH in Norman): {len(t2_predictable)}/{len(t2)}")
print(f"First 20 predictable: {t2_predictable[['gene1','gene2']].head(20).to_dict('records')}")

# Saunders overlap with Norman
sm = open('/tmp/saunders_target_genes.txt').read().splitlines() if Path('/tmp/saunders_target_genes.txt').exists() else []
if sm:
    orth = pd.read_csv('data/external/orthologs/mouse_human_orthologs_symbols.tsv.gz', sep='\t').dropna(subset=['mouse_gene_symbol','human_gene_symbol'])
    orth_o2o = orth[orth['ortholog_type']=='ortholog_one2one']
    m2h = dict(zip(orth_o2o['mouse_gene_symbol'], orth_o2o['human_gene_symbol']))
    saunders_human = sorted({m2h.get(m) for m in sm if m2h.get(m)})
    overlap_norman_saunders = sorted(set(saunders_human) & pert_genes)
    print(f"\n=== NORMAN ∩ SAUNDERS === {len(overlap_norman_saunders)}/{len(saunders_human)} Saunders human-orth genes are in Norman")
    print(f"Shared: {overlap_norman_saunders}")

# Save pert names list for downstream scripts
out_path = Path("Analysis/Perturbation/results/d3_synergy/norman_pert_names.json")
out_path.parent.mkdir(parents=True, exist_ok=True)
with open(out_path, 'w') as fh:
    json.dump(sorted(pert_genes), fh)
print(f"\nWrote {out_path}")

# Final gate decision
if len(overlap_t1) >= 3:
    print(f"\n*** GATE: GREEN ({len(overlap_t1)}/4 tier1 in Norman vocab — proceed with Phase B) ***")
elif len(overlap_t1) >= 2:
    print(f"\n*** GATE: YELLOW ({len(overlap_t1)}/4 tier1 in Norman vocab — borderline, recommend proceed) ***")
else:
    print(f"\n*** GATE: RED ({len(overlap_t1)}/4 tier1 in Norman vocab — HALT, deeper vocab mismatch) ***")

print("\nPRE-FLIGHT DONE")
PY
