"""Recompute each frozen program separately from the pilot integer counts."""
import hashlib
import json
import os
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem
from scipy import sparse

assert os.environ.get('SLURM_JOB_ID'), 'Compute allocation required'
root = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
task = root / 'Analysis/Hypotheses/NUCLEUS-FIX-20261107'
pilot = task / 'pilot_SRX22021596'
for name, expected in json.loads((task / 'checks/preserved_files.json').read_text()).items():
    assert hashlib.sha256((root / name).read_bytes()).hexdigest() == expected, name
membership = pd.read_csv(root / 'Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot/program_membership_v2.tsv', sep='\t')
with h5py.File(pilot / 'unique_counts.h5ad', 'r') as handle:
    obs = read_elem(handle['obs'])
    var = read_elem(handle['var'])
    group = handle['raw/X']
    counts = sparse.csr_matrix((group['data'][:], group['indices'][:], group['indptr'][:]), shape=tuple(group.attrs['shape']))
depth = np.asarray(counts.sum(axis=1)).ravel()
assert (depth > 0).all()
saved = pd.read_csv(pilot / 'fixed_weight_projection_by_capture.tsv.gz', sep='\t')
assert saved.program_uid.nunique() == membership.program_uid.nunique() == 117
errors = []
for uid, members in membership.groupby('program_uid', sort=True):
    j = var.index.get_indexer(members.source_gene)
    assert (j >= 0).all()
    weights = members.original_l1_weight.to_numpy()
    assert abs(weights.sum() - 1) < 1e-12
    # Dense per-program calculation, separate from the producer's sparse
    # 117-column multiplication. Each cell uses its full-transcriptome depth.
    expression = np.log1p(10000 * counts[:, j].toarray() / depth[:, None])
    cell_scores = np.sum(expression * weights[None, :], axis=1)
    for lineage, positions in obs.groupby('cell_type', observed=True).indices.items():
        row = saved[saved.program_uid.eq(uid) & saved.lineage.eq(str(lineage))]
        assert len(row) == 1 and int(row.iloc[0].n_cells) == len(positions)
        error = abs(cell_scores[positions].mean() - float(row.iloc[0].score))
        assert error < 1e-12, (uid, lineage, error)
        errors.append(error)
result = dict(programs=117, cells=len(obs), lineage_program_means_checked=len(errors),
              maximum_absolute_error=max(errors), frozen_files_unchanged=True,
              interpretation='Fixed-weight log1p(CP10k) projections; not native Hotspot graph/PCA scores')
(task / 'checks/pilot_projection_independent.json').write_text(json.dumps(result, indent=2))
print(json.dumps(result, indent=2))
