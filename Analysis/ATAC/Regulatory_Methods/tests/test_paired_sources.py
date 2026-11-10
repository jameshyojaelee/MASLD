"""Independent invariants for source qualification and exact paired RNA joining."""
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / 'scripts/qualify_sources.py'
SPEC = importlib.util.spec_from_file_location('qualify_sources', SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class SourceQualificationTests(unittest.TestCase):
    def test_portal_links_do_not_establish_donor_independence(self):
        result = MODULE.qualification_from_documents({
            'FNIH_readme': {'text': '86 livers https://epigenome.wustl.edu/MASLD/ Raw and supplementary data will be available soon.'},
            'FNIH_portal': {'text': '<a href="https://example.org/matrix.h5ad">counts</a>'}})
        self.assertTrue(result['FNIH']['source_reported_86_livers'])
        self.assertTrue(result['FNIH']['portal_candidate_data_links'])
        self.assertFalse(result['FNIH']['qualified_for_MASH_replication'])
        self.assertEqual(result['FNIH']['independence_from_GSE244832'], 'unverified')

    @unittest.skipUnless(os.environ.get('SLURM_JOB_ID'), 'RNA execution requires compute allocation')
    def test_raw_rna_preserves_requested_cell_order_and_library_namespace(self):
        import h5py
        import numpy as np
        from scipy import sparse
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / 'raw.h5'
            ids = [f'ENSG{i:011d}' for i in range(36601)]
            symbols = [f'G{i}' for i in range(36601)]
            digest = hashlib.sha256()
            for identity, symbol in zip(ids, symbols):
                digest.update(identity.encode() + b'\0' + symbol.encode() + b'\n')
            with h5py.File(source, 'w') as handle:
                matrix = handle.create_group('matrix')
                matrix.create_dataset('data', data=np.array([3, 7, 2], dtype=np.int32))
                matrix.create_dataset('indices', data=np.array([0, 2, 1], dtype=np.int32))
                matrix.create_dataset('indptr', data=np.array([0, 2, 3], dtype=np.int32))
                matrix.create_dataset('shape', data=np.array([36601, 2]))
                matrix.create_dataset('barcodes', data=np.array(['AAA-1', 'BBB-1'], dtype='S'))
                features = matrix.create_group('features')
                features.create_dataset('id', data=np.array(ids, dtype='S'))
                features.create_dataset('name', data=np.array(symbols, dtype='S'))
                features.create_dataset('feature_type', data=np.array(['Gene Expression'] * 36601, dtype='S'))
            lock = root / 'lock.json'
            lock.write_text(json.dumps({'raw_rna_sources': {'well1': {'path': str(source), 'bytes': source.stat().st_size, 'sha256': MODULE.sha(source)}},
                                       'raw_rna_registry': {'rna_feature_id_and_symbol_sha256': digest.hexdigest()}}))
            output = root / 'paired'; output.mkdir()
            with (output / 'cells.tsv').open('w') as handle:
                handle.write('well_id\traw_barcode\tcell_id\nwell1\tBBB-1\twell1_BBB-1\nwell1\tAAA-1\twell1_AAA-1\n')
            with (output / 'genes.tsv').open('w') as handle:
                handle.write('gene_id\tgene_symbol\n')
                handle.writelines(f'{identity}\t{symbol}\n' for identity, symbol in zip(ids, symbols))
            MODULE.export_raw_rna({'project_root': str(root), 'paired': {'source_lock': str(lock)}}, output)
            counts = sparse.load_npz(output / 'rna_counts.npz')
            self.assertEqual(counts[0, 1], 2)
            self.assertEqual(counts[1, 0], 3)
            self.assertEqual(counts[1, 2], 7)
            self.assertEqual(int(counts.sum()), 12)
            with self.assertRaisesRegex(ValueError, 'overwrite'):
                MODULE.export_raw_rna({'paired': {'source_lock': str(lock)}}, output)


if __name__ == '__main__':
    unittest.main()
