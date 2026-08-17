from __future__ import annotations

import unittest

import numpy as np
import pandas as pd

from masld_cl.control_adapter import fit_and_apply_control_offsets
from masld_cl.contracts import ContractError


class TestControlAdapter(unittest.TestCase):
    def fixture(self):
        rows = []
        values = []
        for donor in ("r1", "r2", "r3"):
            rows.append((donor, "REF", "unsorted", "REF_x_unsorted", True, False, False, True))
            values.append([0.0, float(donor[-1])])
        for dataset, shift in (("Q1", 5.0), ("Q2", -4.0)):
            for donor in ("1", "2", "3"):
                rows.append((f"{dataset}c{donor}", dataset, "unsorted", f"{dataset}_x_unsorted", False, True, True, True))
                values.append([shift, float(donor)])
                rows.append((f"{dataset}d{donor}", dataset, "unsorted", f"{dataset}_x_unsorted", False, True, False, True))
                values.append([shift + 2.0, float(donor)])
        cells = pd.DataFrame(rows, columns=[
            "donor_id", "dataset", "preparation", "technical_batch",
            "strict_reference", "primary_query", "query_control", "analysis_eligible",
        ])
        return np.asarray(values, dtype=np.float32), cells

    def test_reference_and_within_batch_differences_are_preserved(self):
        latent, cells = self.fixture()
        result, offsets = fit_and_apply_control_offsets(
            latent, cells, query_datasets=["Q1", "Q2"]
        )
        reference = cells["strict_reference"].to_numpy()
        np.testing.assert_array_equal(result[reference], latent[reference])
        for batch in ("Q1_x_unsorted", "Q2_x_unsorted"):
            mask = cells["technical_batch"].to_numpy() == batch
            np.testing.assert_allclose(
                result[mask] - result[mask][0], latent[mask] - latent[mask][0],
                rtol=0, atol=1e-6,
            )
        self.assertEqual(len(offsets), 2)
        self.assertTrue(all(value["n_control_donors"] == 3 for value in offsets))

    def test_underpowered_controls_fail_closed(self):
        latent, cells = self.fixture()
        keep = ~cells["donor_id"].isin(["Q1c3"])
        with self.assertRaisesRegex(ContractError, "underpowered"):
            fit_and_apply_control_offsets(
                latent[keep], cells.loc[keep].reset_index(drop=True),
                query_datasets=["Q1", "Q2"],
            )


if __name__ == "__main__":
    unittest.main()
