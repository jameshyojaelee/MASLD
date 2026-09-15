from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

try:
    import torch
except ImportError:  # pragma: no cover - exercised only in torch-free contexts
    torch = None

if torch is not None:
    from scripts.b6k_admission_probe import (
        ProbeError,
        arch_flags,
        compare_fixtures,
        execution_path,
        guarded,
        module_matching,
        parse_arch,
        pin_determinism,
        probe_checkpoint_load,
        run_training,
        tensor_digest,
    )


SEED = 20260825

# The two arch lists this campaign actually has, read from the shipped builds.
CU124_ARCH = ["sm_50", "sm_60", "sm_70", "sm_75", "sm_80", "sm_86", "sm_90"]
CU130_ARCH = ["sm_75", "sm_80", "sm_86", "sm_90", "sm_100", "sm_120"]
L40S = (8, 9)
B6K = (12, 0)


@unittest.skipIf(torch is None, "torch is not importable in this interpreter")
class B6KAdmissionProbeTests(unittest.TestCase):
    def test_module_matching_reproduces_dotted_and_numeric_keys(self) -> None:
        state = {
            "encoder.0.weight": torch.randn(4, 3),
            "encoder.0.bias": torch.randn(4),
            "encoder.2.weight": torch.randn(2, 4),
            "head.weight": torch.randn(1, 2),
        }
        module = module_matching(state)
        self.assertEqual(sorted(module.state_dict()), sorted(state))
        result = module.load_state_dict(state, strict=True)
        self.assertEqual(list(result.missing_keys), [])
        self.assertEqual(list(result.unexpected_keys), [])

    def test_module_matching_handles_a_numeric_leaf_parameter_name(self) -> None:
        # The admitted cobolt checkpoint keys look like "beta.0" and
        # "beta_dataset.1", so the leaf itself is numeric rather than the
        # intermediate module.  That is a different code path from "encoder.0".
        state = {
            "beta.0": torch.randn(3),
            "beta.1": torch.randn(3),
            "beta_dataset.0": torch.randn(2),
        }
        module = module_matching(state)
        self.assertEqual(sorted(module.state_dict()), sorted(state))
        result = module.load_state_dict(state, strict=True)
        self.assertEqual(list(result.missing_keys), [])
        self.assertEqual(list(result.unexpected_keys), [])

    def test_module_matching_mixes_float_parameters_and_integer_buffers(self) -> None:
        # cobolt carries both float32 and int64 entries, so a single checkpoint
        # exercises the parameter path and the buffer path at once.
        state = {
            "beta.0": torch.randn(3),
            "counts.0": torch.tensor([1, 2], dtype=torch.int64),
        }
        module = module_matching(state)
        self.assertEqual(
            [name for name, _ in module.named_parameters()], ["beta.0"]
        )
        self.assertEqual([name for name, _ in module.named_buffers()], ["counts.0"])
        module.load_state_dict(state, strict=True)

    def test_module_matching_registers_integer_entries_as_buffers(self) -> None:
        state = {"steps": torch.tensor([3, 4], dtype=torch.int64)}
        module = module_matching(state)
        self.assertEqual([name for name, _ in module.named_buffers()], ["steps"])
        self.assertEqual([name for name, _ in module.named_parameters()], [])

    def test_checkpoint_load_rejects_a_digest_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "state_dict.pt"
            torch.save({"layer.weight": torch.randn(3, 3)}, path)
            with self.assertRaises(ProbeError):
                probe_checkpoint_load(path, "0" * 64)

    def test_checkpoint_load_passes_strictly_on_a_tensor_state_dict(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "state_dict.pt"
            torch.save(
                {"block.0.weight": torch.randn(5, 5), "block.0.bias": torch.randn(5)},
                path,
            )
            record = probe_checkpoint_load(path, None)
        self.assertEqual(record["status"], "pass")
        self.assertEqual(record["tensor_count"], 2)
        self.assertEqual(record["parameter_count"], 30)
        self.assertEqual(record["missing_keys"], [])
        self.assertEqual(record["unexpected_keys"], [])
        self.assertFalse(record["checkpoint_sha256_verified"])

    def test_checkpoint_load_refuses_a_non_tensor_payload(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "state_dict.pt"
            torch.save({"weight": torch.randn(2, 2), "epoch": 7}, path)
            with self.assertRaises(ProbeError):
                probe_checkpoint_load(path, None)

    def test_identical_fixtures_compare_bitwise_identical(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            tensors = {"matmul_fp32": torch.randn(8, 8), "conv2d_fp32": torch.randn(4)}
            left = Path(directory) / "left.pt"
            right = Path(directory) / "right.pt"
            torch.save(tensors, left)
            torch.save(tensors, right)
            report = compare_fixtures(left, right)
        self.assertEqual(report["max_abs_difference"], 0.0)
        self.assertEqual(report["max_rel_difference"], 0.0)
        self.assertTrue(report["all_bitwise_identical"])
        self.assertEqual(report["tensors_compared"], 2)

    def test_comparison_recovers_a_known_perturbation(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            base = torch.ones(4, 4)
            moved = base.clone()
            moved[0, 0] = 1.25
            left = Path(directory) / "left.pt"
            right = Path(directory) / "right.pt"
            torch.save({"matmul_fp32": base}, left)
            torch.save({"matmul_fp32": moved}, right)
            report = compare_fixtures(left, right)
        self.assertAlmostEqual(report["max_abs_difference"], 0.25, places=12)
        self.assertAlmostEqual(report["max_rel_difference"], 0.2, places=12)
        self.assertFalse(report["all_bitwise_identical"])

    def test_comparison_reports_disjoint_tensor_names(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            left = Path(directory) / "left.pt"
            right = Path(directory) / "right.pt"
            torch.save({"shared": torch.zeros(2), "only_left": torch.zeros(2)}, left)
            torch.save({"shared": torch.zeros(2), "only_right": torch.zeros(2)}, right)
            report = compare_fixtures(left, right)
        self.assertEqual(report["left_only"], ["only_left"])
        self.assertEqual(report["right_only"], ["only_right"])
        self.assertEqual(report["tensors_compared"], 1)

    def test_tensor_digest_is_dtype_independent_for_equal_values(self) -> None:
        values = [1.0, 2.5, -3.25]
        self.assertEqual(
            tensor_digest(torch.tensor(values, dtype=torch.float32)),
            tensor_digest(torch.tensor(values, dtype=torch.float64)),
        )

    def test_determinism_knobs_are_all_recorded(self) -> None:
        applied = pin_determinism(SEED)
        expected = {
            "use_deterministic_algorithms",
            "cudnn_deterministic",
            "cudnn_benchmark_off",
            "matmul_tf32_off",
            "cudnn_tf32_off",
            "fp32_matmul_precision_highest",
        }
        self.assertEqual(set(applied), expected)
        # A knob that silently vanished across a torch upgrade would let the
        # cells differ for free, so absence must be recorded, never swallowed.
        for label, outcome in applied.items():
            self.assertTrue(
                outcome == "applied" or outcome.startswith("unavailable:"),
                f"{label} recorded an uninterpretable outcome: {outcome}",
            )

    def test_arch_strings_split_major_and_minor_on_the_last_digit(self) -> None:
        self.assertEqual(parse_arch("sm_86"), (8, 6))
        self.assertEqual(parse_arch("sm_90"), (9, 0))
        # Blackwell is the case a fixed-width slice gets wrong: sm_120 is
        # major 12 minor 0, not major 1 minor 20.
        self.assertEqual(parse_arch("sm_120"), (12, 0))
        self.assertEqual(parse_arch("sm_100"), (10, 0))
        self.assertEqual(parse_arch("sm_120a"), (12, 0))
        self.assertIsNone(parse_arch("compute_90"))
        self.assertIsNone(parse_arch("sm_"))

    def test_b6k_has_no_sass_in_the_pinned_cu124_build(self) -> None:
        # The whole reason the probe campaign exists.
        record = execution_path(B6K, CU124_ARCH)
        self.assertEqual(record["execution_path"], "no_matching_sass")
        self.assertIn("PTX", record["execution_note"])

    def test_b6k_is_native_on_the_cu130_build(self) -> None:
        record = execution_path(B6K, CU130_ARCH)
        self.assertEqual(record["execution_path"], "native_sass")

    def test_l40s_is_recorded_as_sm86_compatibility_not_native(self) -> None:
        # The note the artifact must carry so a later reader does not mistake
        # the L40S cells for native sm_89 execution.
        for arch in (CU124_ARCH, CU130_ARCH):
            record = execution_path(L40S, arch)
            self.assertEqual(
                record["execution_path"], "minor_version_binary_compatibility"
            )
            self.assertEqual(record["execution_served_by"], "sm_86")
            self.assertIn("NOT native sm_89", record["execution_note"])

    def test_arch_flags_reports_a_list_of_targets(self) -> None:
        flags = arch_flags()
        self.assertIsInstance(flags, list)
        self.assertTrue(all(isinstance(flag, str) for flag in flags))

    def test_guarded_converts_a_failure_into_a_recorded_result(self) -> None:
        def explode() -> dict[str, object]:
            raise RuntimeError("no kernel image is available for execution on the device")

        record = guarded(explode)
        self.assertEqual(record["status"], "error")
        self.assertEqual(record["error_type"], "RuntimeError")
        self.assertIn("no kernel image", record["error"])

    def test_training_is_reproducible_for_a_fixed_seed(self) -> None:
        # Probe 5 asserts an exact resume.  That claim is only meaningful if the
        # trajectory is reproducible in the first place, which is checkable on
        # CPU without a GPU node.
        pin_determinism(SEED)
        first, _, _ = run_training(torch.device("cpu"), SEED, 6)
        pin_determinism(SEED)
        second, _, _ = run_training(torch.device("cpu"), SEED, 6)
        self.assertEqual(first, second)
        self.assertLess(first[-1], first[0])


if __name__ == "__main__":
    unittest.main()
