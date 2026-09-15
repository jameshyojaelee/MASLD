from __future__ import annotations

import unittest

from scripts.check_contract_suite_baseline import (
    PREEXISTING_FAILURES,
    ContractBaselineError,
    check,
    parse_failures,
    parse_passed,
    parse_ran,
)


def report(*, failures: tuple[str, ...] = (), ran: int = 80, ok: tuple[str, ...] = ()) -> str:
    # Real unittest -v shape: "method (dotted.path.Class.method) ... ok".
    lines = [f"{name.rsplit(chr(46), 1)[-1]} ({name}) ... ok" for name in ok]
    for name in failures:
        lines.append(
            f"FAIL: {name} (tests.contract.test_registry_files."
            f"RegistryFileContractTest.{name})"
        )
    lines.append(f"Ran {ran} tests in 1.0s")
    return "\n".join(lines)


# Derived from the register, never copied. A hardcoded second copy drifts the
# moment another lane resolves an entry, which is exactly what happened on
# 2026-08-25 when test_model_census_maps_exactly_to_audit_bundles was fixed.
PREEXISTING = tuple(sorted(PREEXISTING_FAILURES)[:1])


class ParseTests(unittest.TestCase):
    def test_failures_and_errors_are_both_collected(self) -> None:
        text = "FAIL: test_a (x.y)\nERROR: test_b (x.y)\nRan 2 tests in 0.1s"
        self.assertEqual(parse_failures(text), {"test_a", "test_b"})

    def test_missing_summary_line_is_an_error(self) -> None:
        with self.assertRaises(ContractBaselineError):
            parse_ran("no summary here")

    def test_passed_names_carry_the_module_path_not_just_the_method(self) -> None:
        text = "test_x (tests.unit.test_mine.Case.test_x) ... ok\nRan 1 tests in 0.1s"
        self.assertEqual(parse_passed(text), {"tests.unit.test_mine.Case.test_x"})


class BaselineCheckTests(unittest.TestCase):
    def test_only_preexisting_failures_pass_the_gate(self) -> None:
        text = report(failures=PREEXISTING, ok=("tests.unit.test_mine.C.test_a",))
        summary = check(text, ["test_mine"], minimum_tests=10)
        self.assertEqual(summary["new_failures"], [])
        self.assertEqual(
            summary["preexisting_failures_still_present"], list(PREEXISTING)
        )

    def test_a_new_failure_fails_the_gate(self) -> None:
        text = report(
            failures=PREEXISTING + ("test_something_i_broke",),
            ok=("tests.unit.test_mine.C.test_a",),
        )
        with self.assertRaises(ContractBaselineError) as caught:
            check(text, ["test_mine"], minimum_tests=10)
        self.assertIn("test_something_i_broke", str(caught.exception))

    def test_a_fixed_preexisting_failure_still_passes(self) -> None:
        """A register entry resolved by another lane must not fail the gate."""

        text = report(failures=(), ok=("tests.unit.test_mine.C.test_a",))
        summary = check(text, ["test_mine"], minimum_tests=10)
        self.assertEqual(
            summary["preexisting_failures_since_fixed"],
            sorted(PREEXISTING_FAILURES),
        )
        self.assertEqual(summary["new_failures"], [])

    def test_a_truncated_run_fails_rather_than_passing_quietly(self) -> None:
        text = report(ran=3, ok=("tests.unit.test_mine.C.test_a",))
        with self.assertRaises(ContractBaselineError):
            check(text, ["test_mine"], minimum_tests=70)

    def test_required_module_must_actually_have_run(self) -> None:
        text = report(ok=("tests.unit.test_other.C.test_a",))
        with self.assertRaises(ContractBaselineError) as caught:
            check(text, ["test_mine"], minimum_tests=10)
        self.assertIn("test_mine", str(caught.exception))


if __name__ == "__main__":
    unittest.main()
