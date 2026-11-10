"""Check called-cell UMI totals against jointly deduplicated molecule records."""
import json
import os

import h5py
import numpy as np
import pandas as pd


def check_molecule_accounting(directory, barcodes, gene_ids, matrix):
    assert os.environ.get('SLURM_JOB_ID')
    by_cell = np.zeros(len(barcodes), dtype=np.int64)
    by_gene = np.zeros(len(gene_ids), dtype=np.int64)
    with h5py.File(directory/'count/outs/molecule_info.h5', 'r') as handle:
        assert np.array_equal(handle['features/id'].asstr()[:], np.asarray(gene_ids))
        assert np.all(np.char.endswith(np.asarray(barcodes, dtype=str), '-1'))
        raw = handle['barcodes'].asstr()[:]
        indices = pd.Index(raw).get_indexer([barcode[:-2] for barcode in barcodes])
        assert (indices >= 0).all() and len(set(indices)) == len(indices)
        lookup = np.full(len(raw), -1, dtype=np.int32)
        lookup[indices] = np.arange(len(indices))
        n_rows = len(handle['barcode_idx'])
        for start in range(0, n_rows, 500000):
            end = min(start+500000, n_rows)
            assert np.all(handle['gem_group'][start:end] == 1)
            assert np.all(handle['library_idx'][start:end] == 0)
            cells = lookup[handle['barcode_idx'][start:end]]
            keep = (cells >= 0) & (handle['count'][start:end] > 0)
            features = handle['feature_idx'][start:end][keep]
            by_cell += np.bincount(cells[keep], minlength=len(by_cell))
            by_gene += np.bincount(features, minlength=len(by_gene))
    assert np.array_equal(by_cell, np.asarray(matrix.sum(axis=1)).ravel())
    assert np.array_equal(by_gene, np.asarray(matrix.sum(axis=0)).ravel())
    result = dict(status='passed', n_molecule_records_checked=n_rows,
                  n_called_barcodes=len(barcodes), n_genes=len(gene_ids), called_cell_umis=int(by_cell.sum()),
                  per_cell_and_per_gene_UMI_counts_exact=True,
                  interpretation='one_positive_read_molecule_record_per_jointly_deduplicated_UMI; include_introns_true')
    output = directory/'molecule_accounting.json'
    if output.exists():
        assert json.loads(output.read_text()) == result
    else:
        output.write_text(json.dumps(result, indent=2))
    return result
