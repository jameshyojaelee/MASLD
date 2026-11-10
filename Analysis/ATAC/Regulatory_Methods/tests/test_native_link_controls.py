import importlib.util
from pathlib import Path
import unittest
import numpy as np
import pandas as pd
spec=importlib.util.spec_from_file_location('controls',Path(__file__).resolve().parents[1]/'scripts/native_link_controls.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

class ReferenceControls(unittest.TestCase):
 def test_weighting_constants_and_whole_RNA_shuffle(self):
  donors=['a','a','b','b','b','b'];rna=pd.DataFrame({'target':[1,2,3,4,5,6],'tf':[4,3,5,6,7,8],'constant':[1]*6})
  motifs=pd.DataFrame({'peak_id':['p','p','p','q'],'tf_gene_id':['tf','tf','target','constant']})
  peaks=pd.DataFrame({'peak_id':['p','q'],'start0':[0,20],'end':[10,30]})
  a=m.tripod_controls(rna,motifs,peaks,0,'target',donors)
  b=m.tripod_controls(rna.iloc[[1,0,5,3,2,4]].reset_index(drop=True),motifs,peaks,0,'target',donors)
  np.testing.assert_allclose(a[['RNA_only','prior_only']],b[['RNA_only','prior_only']],atol=1e-12)
  self.assertEqual(a.n_mapped_TFs.tolist(),[1,1]);self.assertEqual(a.RNA_only.iloc[1],0)
  weights=np.array([.25,.25,.125,.125,.125,.125]);x=rna.tf.to_numpy();y=rna.target.to_numpy();xc=x-weights@x;yc=y-weights@y
  self.assertAlmostEqual(a.RNA_only.iloc[0],abs(np.sum(weights*xc*yc)/np.sqrt(np.sum(weights*xc*xc)*np.sum(weights*yc*yc))))
  altered=m.tripod_controls(rna.assign(target=1,tf=999),motifs,peaks,0,'target',donors)
  np.testing.assert_array_equal(a.prior_only,altered.prior_only)

 def test_bulk_offset_repeat_mapping_and_prior_independence(self):
  identity=pd.DataFrame({'peak_id':['p','p'],'RE_input_position0':[0,1],'prior_column0':[638,640],'expected_identity':['r1','r2'],'actual_identity':['r1','r2'],'match':[True,True],'bulk_TFs':[638,638]})
  layer=np.zeros((2,1180));layer[:,638]=[1,-2];layer[:,640]=[-3,4]
  tf=np.array([[1,8],[2,7],[3,6],[4,5]])
  a=m.linger_controls('g','self',[1,2,3,4],tf,['tf','self'],np.array([[.2,.4],[99,99]]),['a','a','b','b'],identity,layer)
  self.assertAlmostEqual(a.RNA_only.iloc[0],.4);self.assertEqual(a.prior_only.iloc[0],7)
  b=m.linger_controls('g','self',[9,1,3,7],tf*99,['tf','self'],np.array([[.2,.4],[99,99]]),['a','a','b','b'],identity,layer)
  self.assertEqual(a.prior_only.iloc[0],b.prior_only.iloc[0])
  order=[1,0,3,2]
  shuffled=m.linger_controls('g','self',np.array([1,2,3,4])[order],tf[order],['tf','self'],np.array([[.2,.4],[99,99]]),['a','a','b','b'],identity,layer)
  np.testing.assert_allclose(a[['RNA_only','prior_only']],shuffled[['RNA_only','prior_only']],atol=1e-12)
  constant=m.linger_controls('g','self',[1]*4,tf,['tf','self'],np.array([[.2,.4],[99,99]]),['a','a','b','b'],identity,layer)
  self.assertEqual(constant.RNA_only.iloc[0],0)
  bad=identity.assign(bulk_TFs=486)
  with self.assertRaisesRegex(ValueError,'bulkTF'):m.linger_controls('g','self',[1,2,3,4],tf,['tf','self'],np.ones((2,2)),['a','a','b','b'],bad,layer)

if __name__=='__main__':unittest.main()
