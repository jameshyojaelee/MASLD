"""Focused checks of stable molecule sampling and complete record retention."""
import unittest

from audit_native_library_subsample import selected


class MoleculeSamplingTests(unittest.TestCase):
    def test_selection_is_stable_and_rg_specific(self):
        calls = [selected(f'read{i}', 'groupA') for i in range(10000)]
        self.assertEqual(calls, [selected(f'read{i}', 'groupA') for i in range(10000)])
        self.assertNotEqual(calls, [selected(f'read{i}', 'groupB') for i in range(10000)])
        self.assertTrue(any(calls))

    def test_all_mates_hits_and_supplementary_records_share_selection(self):
        name = next(f'read{i}' for i in range(10000) if selected(f'read{i}', 'RG'))
        records = [(name, 'RG', mate, hit, supplementary)
                   for mate in [1, 2] for hit in [1, 2, 3] for supplementary in [False, True]]
        retained = [r for r in records if selected(r[0], r[1])]
        self.assertEqual(retained, records)


if __name__ == '__main__':
    unittest.main()
