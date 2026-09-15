"""Stage 0c's freeze must seal before Arm B is touched, and must disclose Arm A.

Three properties are pinned. The freeze must not compute the correlation
between the two axes -- that restraint is the only thing that makes sealing
first worth anything, and it is easy to break by adding one convenience field.
The Arm A observations this stage descends from must be present by value, or the
derived-after exposure is left for a reader to reconstruct. And the outcome map
must keep the underpowered / not-independent distinction in the outcome name,
because that is the distinction Stage 0b showed decides the interpretation.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
MODULE = ROOT / "scripts" / "freeze_arm_b_axis_count_prespecification.py"
ENDPOINTS = (
    ROOT / "executions" / "model-data-880-21130257-gse135251-source"
    / "source" / "outcomes" / "participant_endpoints.tsv"
)


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


freezer = _load(MODULE, "arm_b_freezer_tested")


class RestraintTests(unittest.TestCase):
    def test_the_freezer_calls_no_correlation_function(self) -> None:
        """Scan the code, not the prose.

        The method block legitimately *describes* partial Spearman, so a
        substring scan over the file matches documentation and proves nothing.
        What matters is that nothing is called. This walks the AST and checks
        every call target.
        """

        import ast

        tree = ast.parse(MODULE.read_text(encoding="utf-8"))
        called: set[str] = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Call):
                target = node.func
                if isinstance(target, ast.Name):
                    called.add(target.id)
                elif isinstance(target, ast.Attribute):
                    called.add(target.attr)
        for forbidden in (
            "spearmanr", "pearsonr", "corrcoef", "rankdata", "cov",
            "spearman_correlation", "unit_ranks", "lstsq", "residualize",
        ):
            self.assertNotIn(forbidden, called, f"the freezer calls {forbidden}")

    def test_the_freezer_imports_no_correlation_machinery(self) -> None:
        import ast

        tree = ast.parse(MODULE.read_text(encoding="utf-8"))
        imported: set[str] = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom):
                imported.update(alias.name for alias in node.names)
            elif isinstance(node, ast.Import):
                imported.update(alias.name for alias in node.names)
        for forbidden in ("rankdata", "spearmanr", "masld_bench.evaluators.metrics"):
            self.assertNotIn(forbidden, imported)

    def test_the_document_carries_no_cross_axis_number(self) -> None:
        """No field may hold a correlation between the two Arm B axes."""

        payload = freezer.build(freezer.read_table(ENDPOINTS))
        blob = repr(payload)
        self.assertNotIn("nas_score|fibrosis_stage", blob)
        self.assertNotIn("fibrosis_stage|nas_score", blob)

    def test_the_restraint_is_asserted_in_the_document(self) -> None:
        payload = freezer.build(freezer.read_table(ENDPOINTS))
        inspected = payload["what_was_inspected_before_freezing"]
        self.assertTrue(inspected["no_gene_expression_was_opened"])
        self.assertTrue(inspected["no_correlation_between_the_two_axes_was_computed"])
        self.assertEqual(sorted(inspected["marginals_and_tie_structure_only"]),
                         ["fibrosis_stage", "nas_score"])


class DisclosureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = freezer.build(freezer.read_table(ENDPOINTS))

    def test_every_parent_is_cited_by_digest(self) -> None:
        parents = self.payload["honest_provenance"]["parents"]
        self.assertEqual(parents, freezer.PARENT_DIGESTS)
        for digest in parents.values():
            self.assertEqual(len(digest), 64)

    def test_arm_a_observations_are_recorded_by_value(self) -> None:
        prior = self.payload["arm_a_prior_observations"]
        self.assertEqual(
            prior["label_spearman_sealed_before_expression_was_opened"][
                "saf_activity_sum|fibrosis"], 0.6494)
        self.assertEqual(
            prior["partial_families_bh_0_05_counts"]["nas_activity_sum|fibrosis"], 805)
        self.assertEqual(
            prior["marginal_bh_0_05_counts"]["fibrosis"], 8)
        self.assertEqual(
            prior["split_half_gene_ordering_reliability_full_length"]["fibrosis"],
            0.3574)

    def test_the_in_cohort_status_is_stated_not_implied(self) -> None:
        self.assertIn(
            "in-cohort",
            self.payload["honest_provenance"]["in_cohort_and_why_that_is_acceptable"],
        )

    def test_the_scale_hazard_is_named_and_not_guessed(self) -> None:
        hazard = self.payload["fibrosis_scale_hazard"]
        self.assertIn("no stage 4", hazard["arm_a_native_scale"])
        self.assertIn(
            "not determinable",
            hazard["whether_arm_a_lacks_a_level_4_or_sampled_no_f4_participant"],
        )
        self.assertIn("no_recode_is_applied_here", hazard)


class OutcomeMapTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = freezer.build(freezer.read_table(ENDPOINTS))

    def test_the_floor_distinction_survives_into_the_outcome_names(self) -> None:
        outcomes = [c["outcome"] for c in self.payload["decision_rule"]["cells"]]
        self.assertIn("ONE_AXIS_ACTIVITY_FIBROSIS_UNDERPOWERED", outcomes)
        self.assertIn("ONE_AXIS_ACTIVITY_FIBROSIS_NOT_INDEPENDENT", outcomes)

    def test_the_map_is_symmetric_over_the_two_directions(self) -> None:
        outcomes = [c["outcome"] for c in self.payload["decision_rule"]["cells"]]
        for suffix in ("UNDERPOWERED", "NOT_INDEPENDENT"):
            self.assertIn(f"ONE_AXIS_ACTIVITY_FIBROSIS_{suffix}", outcomes)
            self.assertIn(f"ONE_AXIS_FIBROSIS_ACTIVITY_{suffix}", outcomes)

    def test_neither_surviving_is_a_declared_cell(self) -> None:
        self.assertIn(
            "NO_AXIS_RESOLVED",
            [c["outcome"] for c in self.payload["decision_rule"]["cells"]],
        )

    def test_e1_is_decisive_and_e2_never_overrides_it(self) -> None:
        criteria = self.payload["criteria"]
        self.assertTrue(criteria["e1_activity_and_fibrosis_are_independent"]["decisive"])
        self.assertFalse(
            criteria["e2_an_empty_direction_is_classified_against_the_measured_floor"][
                "decisive"])
        self.assertTrue(
            self.payload["decision_rule"]["e2_qualifies_the_outcome_name_it_never_overrides_e1"])

    def test_e2_is_a_disambiguation_rule_not_a_pass_fail_gate(self) -> None:
        e2 = self.payload["criteria"][
            "e2_an_empty_direction_is_classified_against_the_measured_floor"]
        self.assertTrue(e2["is_a_disambiguation_rule_not_a_pass_fail_gate"])
        self.assertIn("NO_APPLICABLE_CONDITIONS", e2["applies_only_to_empty_directions"])


class SubstrateTests(unittest.TestCase):
    def test_arm_b_realises_fibrosis_stage_4(self) -> None:
        payload = freezer.build(freezer.read_table(ENDPOINTS))
        counts = payload["what_was_inspected_before_freezing"][
            "marginals_and_tie_structure_only"]["fibrosis_stage"]["value_counts"]
        self.assertEqual(counts, {"0": 35, "1": 41, "2": 48, "3": 44, "4": 12})

    def test_n_is_180_and_not_re_cut(self) -> None:
        payload = freezer.build(freezer.read_table(ENDPOINTS))
        self.assertEqual(payload["cohort"]["n_participants"], 180)
        self.assertIn("216", payload["cohort"]["n_is_180_not_216"])

    def test_a_cohort_without_stage_4_is_refused(self) -> None:
        rows = freezer.read_table(ENDPOINTS)
        capped = [dict(r, fibrosis_stage="3" if r["fibrosis_stage"] == "4" else r["fibrosis_stage"]) for r in rows]
        with self.assertRaises(freezer.PrespecificationError):
            freezer.build(capped)


if __name__ == "__main__":
    unittest.main()
