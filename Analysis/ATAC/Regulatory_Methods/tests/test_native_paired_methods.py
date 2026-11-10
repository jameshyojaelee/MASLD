"""Scientific invariants for bounded native training preparation."""
import importlib.util
from pathlib import Path
import sys
import unittest
import numpy as np
import pandas as pd
from scipy import sparse

scripts=Path(__file__).resolve().parents[1]/'scripts'
sys.path.insert(0,str(scripts))
spec=importlib.util.spec_from_file_location('native',scripts/'native_paired_methods.py')
native=importlib.util.module_from_spec(spec); spec.loader.exec_module(native)


class NativePreparation(unittest.TestCase):
    def test_prior_feature_identity_not_merely_matrix_shape(self):
        tf={10:'TF_A',20:'TF_B',30:'TF_C'}
        re={100:'chr20:100-200',200:'chr20:300-400'}
        rows=native.validate_linger_prior_features([0,2,4],[10,20,30],[100,200],[10,30],[200],tf,re,5)
        self.assertEqual([r['expected_identity'] for r in rows],['TF_A','TF_C','chr20:300-400'])
        # Same lengths and all indices in bounds, but observed-TF offset2
        # selects the wrong RE identity. Shape-only checks would accept it.
        with self.assertRaises(ValueError):
            native.validate_linger_prior_features([0,2,3],[10,20,30],[100,200],[10,30],[200],tf,re,5)
        # The wrong offset can also select a TF column for a RE input.
        with self.assertRaises(ValueError):
            native.validate_linger_prior_features([0,2,2],[10,20,30],[100,200],[10,30],[100],tf,re,5)

    def test_linger_background_contains_distinct_training_rows_only(self):
        meta=pd.DataFrame({'donor_id':np.repeat([str(i) for i in range(20)],2),
                           'metacell_id':[f'train_{i}' for i in range(40)]})
        np.random.seed(42); paired=native.linger_background_indices(40)
        np.random.seed(42); shuffled=native.linger_background_indices(40)
        np.testing.assert_array_equal(paired,shuffled)
        audit=native.validate_linger_background(meta,paired,['heldout'])
        self.assertEqual(set(audit['row_indices']),set(range(40)))
        with self.assertRaises(ValueError): native.validate_linger_background(meta,[0,0],[])
        with self.assertRaises(ValueError): native.validate_linger_background(meta,[40],[])
        with self.assertRaises(ValueError): native.validate_linger_background(meta,[0],['0'])
        np.random.seed(42); expected=np.random.choice(62,50,replace=False); next_expected=np.random.rand()
        np.random.seed(42); actual=native.linger_background_indices(62); next_actual=np.random.rand()
        np.testing.assert_array_equal(actual,expected)
        self.assertEqual(next_actual,next_expected)

    def test_tripod_CPM_is_linear_and_uses_full_axis(self):
        cpm=native.native_tripod_cpm(sparse.csr_matrix([[1.,9.],[2.,18.]]))
        np.testing.assert_allclose(cpm.toarray(),[[100000.,900000.]]*2)
        with self.assertRaises(ValueError): native.native_tripod_cpm(sparse.csr_matrix([[0.,0.]]))

    def test_linger_per_nucleus_units_and_full_RNA_denominator(self):
        rna=sparse.csr_matrix([[1.,0.,99.],[9.,0.,1.]])
        atac=sparse.csr_matrix([[0.,9.],[3.,0.]])
        assignment=sparse.csr_matrix([[.5,.5]])
        mean_rna,mean_atac=native.linger_mean_inputs(rna,atac,assignment)
        self.assertAlmostEqual(mean_rna[0,0],(np.log1p(100)+np.log1p(9000))/2)
        self.assertAlmostEqual(mean_atac[0,0],np.log(2))
        self.assertAlmostEqual(mean_atac[0,1],np.log(10)/2)
        self.assertNotAlmostEqual(mean_rna[0,0],np.log1p(10000*10/110))

    def test_mixed_lineages_are_rejected(self):
        cells=pd.DataFrame({'cohort':['GSE296875']*2,'lineage':['hepatocyte','stellate']})
        with self.assertRaises(ValueError): native.qualified_context(cells)
        self.assertEqual(native.qualified_context(cells.iloc[:1])['lineage'],'hepatocyte')

    def test_no_donor_crosses_any_training_boundary(self):
        donors=np.repeat([str(i) for i in range(10)],4)
        stages=list(native.stages(donors,42))
        self.assertEqual(len(stages),20)
        for name,indices,excluded in stages:
            self.assertFalse(set(donors[indices])&excluded,name)
            self.assertGreaterEqual(len(set(donors[indices])),3)
        again=list(native.stages(donors,42))
        for first,second in zip(stages,again):
            np.testing.assert_array_equal(first[1],second[1])

    def test_missing_linger_priors_never_reports_ready(self):
        result=native.linger_inventory(None,'chr1')
        self.assertEqual(result['status'],'not_run')
        self.assertFalse(any(x['exists'] for x in result['files']))
        self.assertEqual(result['controls']['rna_only'],'unsupported')


if __name__=='__main__': unittest.main()
