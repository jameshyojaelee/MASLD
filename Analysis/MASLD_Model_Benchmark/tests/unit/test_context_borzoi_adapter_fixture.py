from __future__ import annotations

import copy
import json
from pathlib import Path
import tempfile
import unittest

from masld_bench.adapters.context_borzoi_fixture import (
    AdapterFixtureSpec,
    ContextBorzoiFixtureError,
    RNAContextMask,
    SyntheticContextRow,
    derange_context_assignments,
    export_fixture_bundle,
    run_adapter_arm,
    verify_nonpickle_export,
)
from scripts.build_context_borzoi_adapter_fixture import (
    FixtureBuildError,
    build_fixture,
    validate_contract,
)


ROOT = Path(__file__).resolve().parents[2]
CONTRACT = (
    ROOT
    / "config/artifacts/models/context_borzoi/rna_context_adapter_fixture_20260825.json"
)


def context(*values: float | None) -> RNAContextMask:
    return RNAContextMask(
        values=values,
        observed_mask=tuple(value is not None for value in values),
        states=tuple(
            "observed" if value is not None else "structurally_missing"
            for value in values
        ),
    )


def row(
    row_id: str,
    values: tuple[float | None, ...] = (1.0, 2.0, 3.0, 4.0),
    *,
    outer_partition: str = "fold-0",
    atac: bool = False,
) -> SyntheticContextRow:
    return SyntheticContextRow(
        row_id=row_id,
        outer_partition=outer_partition,
        dataset_id="synthetic-dataset",
        assay="synthetic-rna",
        lineage="synthetic-lineage",
        topology="synthetic-unpaired",
        rna=context(*values),
        observed_atac_input=atac,
    )


FEATURES = (
    (0.25, -0.50, 0.75, 1.00),
    (1.25, 1.50, -1.75, 2.00),
)


class ContextBorzoiAdapterFixtureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
        cls.spec = validate_contract(ROOT, cls.contract)

    def test_component_authority_and_trigger_remain_closed(self) -> None:
        self.assertEqual(self.spec.complementarity_trigger_state, "NOT_EVALUATED")
        self.assertEqual(self.spec.parent_component, "UNSELECTED")
        self.assertFalse(self.spec.model_weights_present)
        self.assertFalse(self.spec.architecture_built)
        self.assertFalse(self.spec.training_authorized)
        self.assertFalse(self.spec.prediction_authorized)
        self.assertFalse(self.spec.gpu_execution_authorized)

    def test_missing_mask_is_explicit_and_never_numeric_zero(self) -> None:
        missing = context(1.0, None, 3.0, 4.0)
        missing.validate(4)
        self.assertIsNone(missing.values[1])
        self.assertFalse(missing.observed_mask[1])
        with self.assertRaisesRegex(ContextBorzoiFixtureError, "never zero"):
            RNAContextMask(
                values=(1.0, 0.0, 3.0, 4.0),
                observed_mask=(True, False, True, True),
                states=("observed", "structurally_missing", "observed", "observed"),
            ).validate(4)

    def test_sequence_only_and_missing_context_are_exact_identity(self) -> None:
        target = row("target")
        sequence_only = run_adapter_arm(
            self.spec, FEATURES, target, arm="sequence_only"
        )
        self.assertEqual(sequence_only.values, FEATURES)
        self.assertFalse(sequence_only.rna_context_values_consumed)
        missing = run_adapter_arm(
            self.spec,
            FEATURES,
            row("missing", (1.0, None, 3.0, 4.0)),
            arm="conditioned",
        )
        self.assertEqual(missing.values, FEATURES)
        self.assertEqual(missing.effective_arm, "sequence_only_missing_context")
        self.assertFalse(missing.rna_context_values_consumed)
        self.assertEqual(missing.context_observed_mask, (True, False, True, True))

    def test_trans_only_is_sequence_invariant(self) -> None:
        target = row("target")
        first = run_adapter_arm(self.spec, FEATURES, target, arm="trans_only")
        changed = tuple(tuple(value + 100.0 for value in values) for values in FEATURES)
        second = run_adapter_arm(self.spec, changed, target, arm="trans_only")
        self.assertEqual(first.values, second.values)
        self.assertFalse(first.sequence_values_consumed)
        self.assertTrue(first.rna_context_values_consumed)

    def test_permuted_context_is_deranged_within_complete_stratum(self) -> None:
        rows = (row("a"), row("b", (4.0, 3.0, 2.0, 1.0)))
        assignments = derange_context_assignments(
            rows, width=4, seed=self.spec.permutation_seed
        )
        for target in rows:
            source = assignments[target.row_id]
            self.assertNotEqual(source.row_id, target.row_id)
            self.assertEqual(source.permutation_stratum(), target.permutation_stratum())
        conditioned = run_adapter_arm(
            self.spec, FEATURES, rows[0], arm="conditioned"
        )
        permuted = run_adapter_arm(
            self.spec,
            FEATURES,
            rows[0],
            arm="permuted_context",
            permuted_source=assignments[rows[0].row_id],
        )
        self.assertNotEqual(conditioned.values, permuted.values)

    def test_observed_atac_stays_in_separate_family_native_task(self) -> None:
        with self.assertRaisesRegex(ContextBorzoiFixtureError, "observed ATAC"):
            run_adapter_arm(
                self.spec, FEATURES, row("atac-row", atac=True), arm="conditioned"
            )

    def test_export_is_json_only_and_contains_no_model_state(self) -> None:
        target = row("target")
        source = row("source", (4.0, 3.0, 2.0, 1.0))
        assignments = {target.row_id: source, source.row_id: target}
        outputs = {
            "conditioned": run_adapter_arm(
                self.spec, FEATURES, target, arm="conditioned"
            ),
            "sequence_only": run_adapter_arm(
                self.spec, FEATURES, target, arm="sequence_only"
            ),
            "trans_only": run_adapter_arm(
                self.spec, FEATURES, target, arm="trans_only"
            ),
            "permuted_context": run_adapter_arm(
                self.spec,
                FEATURES,
                target,
                arm="permuted_context",
                permuted_source=source,
            ),
            "missing_context": run_adapter_arm(
                self.spec,
                FEATURES,
                row("missing", (1.0, None, 3.0, 4.0)),
                arm="conditioned",
            ),
        }
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "export"
            manifest = export_fixture_bundle(output, self.spec, outputs, assignments)
            self.assertFalse(manifest["pickle_present"])
            self.assertFalse(manifest["model_weights_present"])
            self.assertFalse(manifest["tensor_state_present"])
            self.assertEqual(
                {path.name for path in output.iterdir()},
                {
                    "COMPLETE",
                    "export_manifest.json",
                    "fixture_contract.json",
                    "fixture_outputs.json",
                },
            )
            (output / "rogue.pkl").write_bytes(b"\x80\x04fixture")
            with self.assertRaises(ContextBorzoiFixtureError):
                verify_nonpickle_export(output)

    def test_contract_cannot_authorize_architecture_or_gpu(self) -> None:
        changed = copy.deepcopy(self.contract)
        changed["conditional_gate"]["architecture_authorized"] = True
        with self.assertRaises(FixtureBuildError):
            validate_contract(ROOT, changed)
        changed = copy.deepcopy(self.contract)
        changed["execution_disposition"]["gpu_execution_authorized"] = True
        with self.assertRaises(FixtureBuildError):
            validate_contract(ROOT, changed)

    def test_end_to_end_builder_uses_only_synthetic_inputs(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "fixture"
            receipt = build_fixture(ROOT, self.contract, output)
            self.assertEqual(
                receipt["status"], "pass_synthetic_fixture_architecture_blocked"
            )
            self.assertFalse(receipt["biological_data_read"])
            self.assertFalse(receipt["observed_atac_read"])
            self.assertFalse(receipt["model_weights_present"])
            self.assertFalse(receipt["gpu_executed"])
            self.assertFalse(receipt["conditional_model_authorized"])


if __name__ == "__main__":
    unittest.main()
