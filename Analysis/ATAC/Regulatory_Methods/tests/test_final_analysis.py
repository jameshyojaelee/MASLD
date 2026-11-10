#!/usr/bin/env python3
"""Scientific invariants for retained-nucleus controls and training omission."""
import importlib.util
import ast
import os
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd
from scipy import sparse

A=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(A/'scripts'))
import final_prediction_refits as refits
from final_rna_molecule_qc import check_molecule_accounting
import h5py


class ScientificInvariants(unittest.TestCase):
    def test_molecule_accounting_detects_wrong_gene_or_cell_with_same_total(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); (root/'count/outs').mkdir(parents=True)
            with h5py.File(root/'count/outs/molecule_info.h5','w') as handle:
                string=h5py.string_dtype()
                handle.create_dataset('features/id',data=['g1','g2'],dtype=string)
                handle.create_dataset('barcodes',data=['AA','BB','CC'],dtype=string)
                for key,value in dict(barcode_idx=[0,0,1,2],feature_idx=[0,1,0,1],
                                      gem_group=[1]*4,library_idx=[0]*4,count=[3,4,2,5]).items():
                    handle.create_dataset(key,data=value)
            correct=sparse.csr_matrix([[1,1],[1,0]])
            result=check_molecule_accounting(root,['AA-1','BB-1'],['g1','g2'],correct)
            self.assertEqual(result['called_cell_umis'],3)
            # Preserve the grand total but put UMIs in the wrong genes or cells.
            for wrong in [[[0,2],[0,1]],[[0,0],[2,1]]]:
                with self.assertRaises(AssertionError):
                    check_molecule_accounting(root,['AA-1','BB-1'],['g1','g2'],sparse.csr_matrix(wrong))

    def test_cicero_graph_keeps_original_indices_after_zero_filtering(self):
        refits.load_core()
        peaks = pd.DataFrame({'peak_id':['promoter', 'distal']}, index=[2, 5])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'links.tsv'
            pd.DataFrame(dict(peak1=['promoter'], peak2=['distal'], coaccess=[.8])).to_csv(path, sep='\t', index=False)
            graph = refits.CORE.CiceroGraph.from_tsv(path, peaks)
            self.assertEqual(graph.csr.shape, (6, 6))
            self.assertEqual(graph.promoter_scores(np.array([2])), {5:.8})
            self.assertEqual(graph.promoter_scores(np.array([0])), {})

    def test_existing_coordinate_and_fold_fixtures(self):
        # These existing fixtures use only assert/numpy. Execute their original
        # bodies without adding pytest to the analysis environment.
        names = {
            'test_promoters_exclude_target_features_on_both_strands',
            'test_indexed_geometry_matches_exact_interval_masks',
            'test_fold_seed_is_identical_for_python_and_numpy_scalar_ids',
        }
        tree = ast.parse((A/'tests/test_network.py').read_text())
        functions = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name in names]
        self.assertEqual(len(functions), len(names))
        refits.load_core()
        namespace = dict(module=refits.CORE, np=np, pd=pd)
        exec(compile(ast.Module(body=functions, type_ignores=[]), 'existing_scientific_fixtures', 'exec'), namespace)
        for name in sorted(names):
            with self.subTest(fixture=name):
                namespace[name]()

    def test_shuffle_retains_all_molecules_in_evaluated_donors(self):
        cells=pd.DataFrame(dict(cohort=['c']*220,lineage=['hepatocyte']*220,donor_id=['a']*110+['b']*110))
        groups=[np.arange(0,50),np.arange(50,100),np.arange(110,160),np.arange(160,210)]
        rna=sparse.csr_matrix(np.arange(660).reshape(220,3))
        retained=np.concatenate(groups)
        for seed in [43,44,45]:
            perm=refits.retained_permutation(cells,groups,seed)
            self.assertEqual(set(perm[retained]),set(retained))
            self.assertTrue(np.array_equal(cells.donor_id.to_numpy(),cells.iloc[perm].donor_id.to_numpy()))
            self.assertTrue(np.array_equal(perm[np.setdiff1d(np.arange(220),retained)],np.setdiff1d(np.arange(220),retained)))
            for donor in ['a','b']:
                idx=retained[cells.iloc[retained].donor_id.eq(donor)]
                self.assertTrue(np.array_equal(rna[idx].sum(axis=0),rna[perm[idx]].sum(axis=0)))
            self.assertFalse(np.array_equal(perm[retained],retained))

    def test_omission_keeps_test_donors_and_inner_assignments(self):
        donors=np.repeat([str(i) for i in range(20)],3)
        outer={str(i):i%5 for i in range(20)}
        inner={str(f):{str(i):i%3 for i in range(20) if i%5!=f} for f in range(5)}
        for f in range(5):
            _,train,test,folds=refits.fixed_donor_splits(donors,outer,inner,f)
            _,reduced,retest,refolds=refits.fixed_donor_splits(donors,outer,inner,f,'7')
            self.assertTrue(np.array_equal(test,retest))
            self.assertNotIn('7',donors[reduced])
            keep=donors[train]!='7'
            self.assertTrue(np.array_equal(train[keep],reduced))
            self.assertTrue(np.array_equal(folds[keep],refolds))
            if f==outer['7']:
                self.assertIn('7',donors[retest])
                self.assertTrue(np.array_equal(train,reduced))


if __name__=='__main__':
    assert os.environ.get('SLURM_JOB_ID'),'Run tests on compute'
    unittest.main()
