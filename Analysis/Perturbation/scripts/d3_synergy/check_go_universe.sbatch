#!/bin/bash
#SBATCH --job-name=GEARS
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/Perturbation/scripts/d3_synergy/logs/check_go_universe_%j.log

set -o pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export PYTHONNOUSERSITE=1
unset PYTHONPATH || true

eval "$(micromamba shell hook --shell bash)"
micromamba activate -p .envs/perturbation_gears

set -eu

python - <<'PY'
import pickle, json
from pathlib import Path
import pandas as pd
import numpy as np

# GO-derived perturbable-gene set (this is what GEARS uses as pert_names base)
p = Path("Analysis/Perturbation/results/finetuned_checkpoints/essential_all_data_pert_genes.pkl")
print(f"Loading {p}, size {p.stat().st_size/1e6:.2f} MB")
with open(p, "rb") as fh:
    data = pickle.load(fh)
print(f"type: {type(data).__name__}; shape/len: {getattr(data,'shape',len(data))}")
if isinstance(data, np.ndarray):
    pert_set = set(data.tolist())
elif isinstance(data, (list, set)):
    pert_set = set(data)
elif isinstance(data, dict):
    pert_set = set(data.keys())
else:
    pert_set = set()
pert_set.discard("ctrl")
print(f"unique non-ctrl entries: {len(pert_set)}")
print(f"sample 20: {sorted(pert_set)[:20]}")

t1 = {"TADA1","VPS41","HECTD1","RNF103"}
print(f"\n=== essential_all_data_pert_genes ({len(pert_set)}) ===")
print(f"Tier1 overlap: {len(t1 & pert_set)}/4 -> {sorted(t1 & pert_set)}")
t2 = pd.read_csv("Analysis/Perturbation/data/hits/d3_tier2_pairs_invocab.csv")
t2g = set(t2.gene1) | set(t2.gene2)
print(f"Tier2 unique genes overlap: {len(t2g & pert_set)}/{len(t2g)}")
both = t2[t2.gene1.isin(pert_set) & t2.gene2.isin(pert_set)]
print(f"Tier2 PAIRS predictable: {len(both)}/{len(t2)}")
print(f"First 10 predictable tier2 pairs: {both.head(10)[['gene1','gene2']].to_dict('records')}")

# Save the pert universe for the union-finetune script
with open("Analysis/Perturbation/results/d3_synergy/gears_essential_pert_universe.json","w") as fh:
    json.dump(sorted(pert_set), fh)
print(f"\nWrote gears_essential_pert_universe.json")

# Saunders panel x GEARS essential_all overlap
sm = open("/tmp/saunders_target_genes.txt").read().splitlines() if Path("/tmp/saunders_target_genes.txt").exists() else []
if sm:
    orth = pd.read_csv("data/external/orthologs/mouse_human_orthologs_symbols.tsv.gz", sep="\t").dropna(subset=["mouse_gene_symbol","human_gene_symbol"])
    orth_o2o = orth[orth["ortholog_type"]=="ortholog_one2one"]
    m2h = dict(zip(orth_o2o["mouse_gene_symbol"], orth_o2o["human_gene_symbol"]))
    saunders_human = sorted({m2h.get(m) for m in sm if m2h.get(m)})
    overlap = set(saunders_human) & pert_set
    print(f"\nSaunders ({len(saunders_human)}) ∩ essential_all_pert_genes ({len(pert_set)}): {len(overlap)}")
    print(f"  shared genes: {sorted(overlap)[:30]}...")
PY
