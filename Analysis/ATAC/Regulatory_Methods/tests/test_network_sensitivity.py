import importlib.util
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

path = Path(__file__).resolve().parents[1] / 'scripts/network_sensitivity.py'
spec = importlib.util.spec_from_file_location('sensitivity', path)
s = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s)


class SensitivityTests(unittest.TestCase):
    def test_pair_shuffle_is_exact_donor_lineage_permutation(self):
        cells = pd.DataFrame({'cohort': ['c'] * 24, 'donor_id': ['a'] * 12 + ['b'] * 12, 'lineage': ['hep'] * 6 + ['mes'] * 6 + ['hep'] * 6 + ['mes'] * 6})
        order = s.shuffle_index(cells, 42)
        self.assertEqual(sorted(order), list(range(24)))
        pd.testing.assert_frame_equal(cells.iloc[order].reset_index(drop=True), cells)
        np.testing.assert_array_equal(order, s.shuffle_index(cells, 42))
        self.assertFalse(np.array_equal(order, np.arange(24)))

    def test_omission_plan_scopes_reused_ids_to_one_cohort(self):
        from unittest.mock import patch
        import json
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary)
            inputs=root/'inputs'; inputs.mkdir()
            pd.DataFrame({'cohort':['a','b'],'donor_id':['1','1'],'lineage':['hep','hep']}).to_csv(inputs/'cells.tsv',sep='\t',index=False)
            baseline=root/'baseline';baseline.mkdir()
            (baseline/'completion.json').write_text('{"status":"completed"}')
            for cohort in ['a','b']:
                context=baseline/(cohort+'_hep');context.mkdir()
                pd.DataFrame({'donor_id':['1'],'outer_fold':[0]}).to_csv(context/'fold_assignments.tsv',sep='\t',index=False)
            (root/'config.json').write_text(json.dumps({'network': {
                'external_stage_links': {'native': {'0:outer': '/old/links.tsv'}},
                'external_stage_roots': {'native': '/old'},
                'external_method_controls': {'native': 'shuffle'}}}))
            (root/'genes.tsv').write_text('gene_id\n')
            arguments=['sensitivity','--config',str(root/'config.json'),'--input-root',str(inputs),'--genes',str(root/'genes.tsv'),'--baseline',str(baseline),'--output',str(root/'out'),'--stages','omission']
            with patch('sys.argv',arguments):s.main()
            plan=json.loads((root/'out/sensitivity_plan.json').read_text())
            self.assertEqual(len(plan['scenarios']),2)
            self.assertEqual({x['omit_cohort'] for x in plan['scenarios']},{'a','b'})
            self.assertEqual(len({x['name'] for x in plan['scenarios']}),2)
            self.assertEqual(plan['excluded_external_methods'], ['native'])

    def test_refit_preserves_outer_donor_folds_and_paired_sums(self):
        from types import SimpleNamespace
        from scipy import sparse
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            baseline = root / 'baseline'
            baseline.mkdir()
            pd.DataFrame({'donor_id':['a','b'], 'outer_fold':[2,4]}).to_csv(baseline/'fold_assignments.tsv',sep='\t',index=False)
            metadata = pd.DataFrame({'donor_id':['a','b']})
            captured = {}
            def build(*args):
                return sparse.eye(2,format='csr'), metadata, [np.array([0,1]), np.array([2,3])], metadata
            def fit(counts, meta, peaks, genes, output, settings, seed, rna, rna_genes):
                output.mkdir()
                captured.update(settings=settings, rna=rna.toarray())
            module = SimpleNamespace(build_metacells=build, fit_group=fit)
            s.refit(module, None, None, None, None, sparse.csr_matrix([[1],[2],[4],[8]]), None, {}, 42, baseline, root/'result')
            self.assertEqual(captured['settings']['outer_fold_map'], {'a':2,'b':4})
            np.testing.assert_array_equal(captured['rna'], [[3],[12]])

    def test_stability_coverage_and_donor_averaging(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            a, b = root / 'a', root / 'b'
            a.mkdir(); b.mkdir()
            pd.DataFrame({'outer_fold':[0,0], 'gene_id':['g','g'], 'method':['m','m'], 'peak_id':['p1','p2'], 'score':[1,2]}).to_csv(a/'links.tsv',sep='\t',index=False)
            pd.DataFrame({'outer_fold':[0,0], 'gene_id':['g','g'], 'method':['m','m'], 'peak_id':['p2','p3'], 'score':[3,4]}).to_csv(b/'links.tsv',sep='\t',index=False)
            common = {'outer_fold':[0,0], 'gene_id':['g','h'], 'endpoint':['rna','rna'], 'method':['m','m'], 'donor_id':['d','d']}
            pd.DataFrame(dict(common,squared_error=[4,8])).to_csv(a/'predictions.tsv',sep='\t',index=False)
            pd.DataFrame(dict(common,squared_error=[2,4])).to_csv(b/'predictions.tsv',sep='\t',index=False)
            s.compare_tables(a,b,root/'out')
            edges=pd.read_csv(root/'out/edge_stability.tsv',sep='\t')
            self.assertAlmostEqual(edges.top_edge_jaccard.iloc[0],1/3)
            donor=pd.read_csv(root/'out/donor_prediction_difference.tsv',sep='\t')
            self.assertEqual(len(donor),1)
            self.assertEqual(donor.baseline_minus_sensitivity_squared_error.iloc[0],3)
            self.assertEqual(donor.n_common_gene_fold_pairs.iloc[0],2)
            for directory in [a,b]:
                table=pd.read_csv(directory/'predictions.tsv',sep='\t')
                extra=table.assign(donor_id='e')
                pd.concat([table,extra]).to_csv(directory/'predictions.tsv',sep='\t',index=False)
            s.compare_tables(a,b,root/'out_two_donors')
            bootstrap=pd.read_csv(root/'out_two_donors/donor_bootstrap_difference.tsv',sep='\t')
            self.assertEqual(bootstrap.n_donors.iloc[0],2)
            self.assertEqual(bootstrap.bootstrap_reps.iloc[0],2000)
            self.assertEqual(bootstrap.ci_low.iloc[0],3)
            self.assertEqual(bootstrap.ci_high.iloc[0],3)

    def test_not_refitted_method_is_distinct_from_refitted_gene_coverage_loss(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            a, b = root/'a', root/'b'
            a.mkdir(); b.mkdir()
            edges = pd.DataFrame({'outer_fold':[0,0,0], 'gene_id':['g','h','g'],
                'method':['m','m','native'], 'peak_id':['p1','p2','p3'], 'score':[1,2,3]})
            edges.to_csv(a/'links.tsv',sep='\t',index=False)
            edges.iloc[[0]].to_csv(b/'links.tsv',sep='\t',index=False)
            pred = pd.DataFrame({'outer_fold':[0,0], 'gene_id':['g','g'], 'endpoint':['rna','rna'],
                'method':['m','native'], 'donor_id':['d','d'], 'squared_error':[1,2]})
            pred.to_csv(a/'predictions.tsv',sep='\t',index=False)
            pred.iloc[[0]].to_csv(b/'predictions.tsv',sep='\t',index=False)
            s.compare_tables(a,b,root/'out')
            result = pd.read_csv(root/'out/edge_stability.tsv',sep='\t')
            excluded = result[result.method == 'native'].iloc[0]
            self.assertEqual(excluded.coverage_state,'method_not_refit')
            self.assertEqual(excluded.baseline_edges,1)
            self.assertEqual(excluded.sensitivity_edges,0)
            self.assertTrue(pd.isna(excluded.top_edge_jaccard))
            self.assertTrue(pd.isna(excluded.common_edge_spearman))
            lost = result[(result.method == 'm') & (result.gene_id == 'h')].iloc[0]
            self.assertEqual(lost.coverage_state,'baseline_only')
            self.assertEqual(lost.top_edge_jaccard,0)



if __name__ == '__main__':
    unittest.main()
