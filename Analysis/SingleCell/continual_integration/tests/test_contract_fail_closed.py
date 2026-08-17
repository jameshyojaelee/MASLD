from __future__ import annotations

import unittest
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from masld_cl.contracts import ContractError, DeterministicGzipTextWriter, sha256_path
from masld_cl.data import compound_cell_ids, verify_prepared


class TestContractFailClosed(unittest.TestCase):
    def test_gzip_manifests_are_byte_deterministic(self):
        with tempfile.TemporaryDirectory() as directory:
            paths = [Path(directory) / f"m{i}.tsv.gz" for i in range(2)]
            for path in paths:
                with DeterministicGzipTextWriter(path) as handle:
                    handle.write("a\tb\n1\t2\n")
            self.assertEqual(sha256_path(paths[0]), sha256_path(paths[1]))

    def test_duplicate_compound_cell_ids_fail(self):
        with self.assertRaisesRegex(ContractError, "not unique"):
            compound_cell_ids(np.array(["L1", "L1"]), np.array(["C1", "C1"]))

    def _fixture(self):
        obs = pd.DataFrame({
            "library_id": ["L1", "L2"], "assay_id": ["A1", "A2"],
            "donor_id": ["D1", "D2"],
            "native_condition_authoritative": ["Healthy", "MASH"],
            "harmonized_stage_authoritative": ["Healthy", "Steatohepatitis"],
            "technical_batch": ["S|p", "S|p"], "audit_cell_type": ["A", "B"],
            "strict_reference": [True, False], "primary_query": [False, True],
            "query_control": [False, False], "analysis_eligible": [True, True],
        }, index=["L1|C1", "L2|C2"])
        matrix = sparse.csr_matrix(np.array([[1, 0], [0, 2]], dtype=np.int32))
        value = ad.AnnData(matrix, obs=obs, var=pd.DataFrame(index=["g1", "g2"]))
        value.layers["counts"] = matrix.copy()
        return value

    def test_noninteger_negative_and_gene_permutation_fail(self):
        value = self._fixture()
        contract = {"observed_contract": {"descriptive": {"cells": 2}}, "raw_counts_sha256": "raw"}
        value.uns["masld_cl_contract"] = {
            "descriptive_cells": 2, "config_sha256": "cfg", "raw_counts_sha256": "raw",
            "model_gene_order_sha256": "wrong", "prepared_obs_sha256": "obs",
            "prepared_counts_sha256": "counts",
        }
        config = {
            "features": {"n_hvg": 2, "selection_batch_key": "technical_batch"},
            "input": {"forbidden_conditioning_fields": ["stage"]}, "_config_sha256": "cfg",
        }
        prepared = {
            "model_gene_order_sha256": "wrong", "prepared_obs_sha256": "obs",
            "prepared_counts_sha256": "counts",
        }
        with self.assertRaisesRegex(ContractError, "genes are missing or permuted"):
            verify_prepared(value, config, contract, prepared)
        # Move past the gene check and exercise count validation.
        from masld_cl.data import _hash_strings
        gene_hash = _hash_strings(value.var_names)
        value.uns["masld_cl_contract"]["model_gene_order_sha256"] = gene_hash
        prepared["model_gene_order_sha256"] = gene_hash
        value.layers["counts"] = sparse.csr_matrix(np.array([[1.5, 0], [0, -2.0]]))
        with self.assertRaisesRegex(ContractError, "noninteger, negative, or non-finite"):
            verify_prepared(value, config, contract, prepared)
