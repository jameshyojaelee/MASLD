import importlib.util
import itertools
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd

SCRIPT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPT_ROOT))
source = Path(os.environ.get('MASLD_EXACT_PERM_SOURCE', SCRIPT_ROOT / '11_run_null_permutations.py'))
spec = importlib.util.spec_from_file_location('exact_permutations', source)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class ExactPermutationTests(unittest.TestCase):
    def test_sign_flip_enumerates_observed_once(self):
        values = np.array([1., 2., 3.])
        null, probability = module.exact_sign_flip(values)
        self.assertEqual(len(null), 8)
        self.assertEqual(probability, 2 / 8)
        oracle = [sum(v * s for v, s in zip(values, signs)) / len(values)
                  for signs in itertools.product([-1, 1], repeat=3)]
        np.testing.assert_array_equal(sorted(null), sorted(oracle))

    def test_ties_and_undefined_inputs(self):
        self.assertEqual(module.exact_sign_flip(np.zeros(3))[1], 1)
        for values in ([], [1, np.nan], [1, np.inf], [[1, 2]]):
            with self.assertRaises(ValueError):
                module.exact_sign_flip(values)

    def test_biopsy_exact_label_oracle(self):
        rows, metadata = [], []
        for participant, delta in enumerate([0., 1., 2., 3.]):
            for timepoint, score in [('baseline', 0.), ('followup', delta)]:
                sample = f'p{participant}_{timepoint}'
                rows.append(dict(sample_id=sample, scoring_scheme='weighted', program_uid='program1', program_score=score))
                metadata.append(dict(sample_id=sample, dataset_id='GSE106737', biological_unit_id=f'p{participant}',
                                     condition='lifestyle_responder' if participant < 2 else 'lifestyle_nonresponder',
                                     timepoint=timepoint))
        with tempfile.TemporaryDirectory(prefix='masld_exact_perm_') as tmp:
            root = Path(tmp)
            path = root / 'analyses/GSE106737/per_sample_program_scores.tsv'
            path.parent.mkdir(parents=True)
            pd.DataFrame(rows).to_csv(path, sep='\t', index=False)
            with patch.object(module, 'CANDIDATE_ROOT', root), patch.object(module, 'read_tsv', return_value=metadata):
                result = module.biopsy({'program1'}, 'fixture')[0]
        self.assertEqual(result['observed'], -2)
        self.assertEqual(result['n_biological'], 4)
        self.assertEqual(result['n_permutations'], 6)
        self.assertEqual(result['empirical_p'], 2 / 6)


if __name__ == '__main__':
    unittest.main()
