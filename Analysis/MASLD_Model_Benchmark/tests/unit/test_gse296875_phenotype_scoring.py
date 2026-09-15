"""Prove the per-fold metric path is unreachable, not merely discouraged."""

from __future__ import annotations

import inspect
from pathlib import Path
import tempfile
import unittest

from masld_bench import gse296875_phenotype_scoring as scoring
from masld_bench.gse296875_phenotype_scoring import (
    EndpointRoster,
    PooledMetricViolation,
    PredictionTableError,
    UnknownScope,
    confirmatory_family,
    endpoint_values,
    random_score_reference,
    load_predictions,
    permutation_null,
    pooled_donor_bootstrap,
    pooled_metric,
    registered_roster,
)


ROOT = Path(__file__).resolve().parents[2]
ENDPOINTS = ROOT / "executions/gse296875-phenotype-endpoints-20260825/endpoint_lock/donor_endpoints.tsv"
FOLDS = ROOT / "executions/gse296875-donor-folds-20260825/split_lock/donor_folds.tsv"

#: The four donors held out in outer fold 3.  Zero of them are fibrosis
#: positive, so a within-fold AUPRC is undefined for this fold.
FOLD_3_DONORS = ("331", "366", "382", "485")


def _rows(path: Path) -> list[dict[str, str]]:
    return scoring.read_tsv(path)


class ClosedScopeRegistryTests(unittest.TestCase):
    def test_no_per_fold_scope_is_registered(self) -> None:
        for name in scoring.SCOPES:
            self.assertNotIn("fold", name)

    def test_fold_scope_cannot_be_requested(self) -> None:
        endpoint_rows = _rows(ENDPOINTS)
        fold_rows = _rows(FOLDS)
        for name in ("outer_fold_3", "fold3", "per_fold", ""):
            with self.assertRaises(UnknownScope):
                registered_roster("fibrosis", name, endpoint_rows, fold_rows)

    def test_registry_is_immutable(self) -> None:
        self.assertIsInstance(scoring.SCOPES, frozenset)
        with self.assertRaises(AttributeError):
            scoring.SCOPES.add("outer_fold_3")  # type: ignore[attr-defined]

    def test_roster_rejects_unregistered_scope_at_construction(self) -> None:
        with self.assertRaises(UnknownScope):
            EndpointRoster(
                endpoint="fibrosis",
                scope="outer_fold_3",
                donor_ids=frozenset(FOLD_3_DONORS),
            )


class PooledOnlyEnforcementTests(unittest.TestCase):
    def setUp(self) -> None:
        self.endpoint_rows = _rows(ENDPOINTS)
        self.fold_rows = _rows(FOLDS)
        self.fibrosis = endpoint_values("fibrosis", self.endpoint_rows)
        self.steatosis = endpoint_values("steatosis", self.endpoint_rows)
        self.roster = registered_roster(
            "fibrosis", "pooled_all_donors", self.endpoint_rows, self.fold_rows
        )
        self.predictions = {
            donor: float(index) for index, donor in enumerate(sorted(self.roster.donor_ids))
        }

    def test_frozen_roster_sizes_match_the_taskspec(self) -> None:
        self.assertEqual(len(self.roster.donor_ids), 37)
        steatosis = registered_roster(
            "steatosis", "pooled_all_donors", self.endpoint_rows, self.fold_rows
        )
        self.assertEqual(len(steatosis.donor_ids), 38)
        self.assertEqual(sum(self.fibrosis.values()), 15.0)
        self.assertEqual(len(self.fibrosis) - sum(self.fibrosis.values()), 22.0)

    def test_single_fold_donors_raise(self) -> None:
        fold_only = {
            donor: self.predictions[donor]
            for donor in FOLD_3_DONORS
            if donor in self.predictions
        }
        self.assertEqual(len(fold_only), 4)
        with self.assertRaises(PooledMetricViolation):
            pooled_metric("auprc", fold_only, self.fibrosis, self.roster)

    def test_fold_3_has_no_fibrosis_positive(self) -> None:
        positives = sum(self.fibrosis[donor] for donor in FOLD_3_DONORS)
        self.assertEqual(positives, 0.0)

    def test_any_proper_subset_raises(self) -> None:
        donors = sorted(self.roster.donor_ids)
        subset = {donor: self.predictions[donor] for donor in donors[:-1]}
        with self.assertRaises(PooledMetricViolation):
            pooled_metric("auprc", subset, self.fibrosis, self.roster)

    def test_masked_donor_superset_raises(self) -> None:
        superset = dict(self.predictions)
        superset["151"] = 0.5
        with self.assertRaises(PooledMetricViolation):
            pooled_metric("auprc", superset, self.fibrosis, self.roster)

    def test_exact_roster_is_accepted(self) -> None:
        value = pooled_metric("auprc", self.predictions, self.fibrosis, self.roster)
        self.assertGreater(value, 0.0)
        self.assertLessEqual(value, 1.0)

    def test_no_per_fold_metric_function_exists(self) -> None:
        source = inspect.getsource(scoring)
        for banned in ("groupby", "by_fold", "per_fold", "fold_metric", "mean_of_folds"):
            self.assertNotIn(banned, source)


class PredictionProjectionTests(unittest.TestCase):
    TABLE = (
        "donor_id\tarm\tlineage_scope\tendpoint\tprediction\touter_fold\n"
        "174\tmolecular\tall_lineage\tfibrosis\t0.10\t0\n"
        "331\tmolecular\tall_lineage\tfibrosis\t0.90\t3\n"
        "174\tmetadata\tall_lineage\tfibrosis\t0.40\t0\n"
    )

    def _write(self, text: str) -> Path:
        handle = tempfile.NamedTemporaryFile(
            "w", suffix=".tsv", delete=False, encoding="utf-8"
        )
        handle.write(text)
        handle.close()
        return Path(handle.name)

    def test_loader_returns_no_fold_information(self) -> None:
        path = self._write(self.TABLE)
        loaded = load_predictions(
            path, arm="molecular", lineage_scope="all_lineage", endpoint="fibrosis"
        )
        self.assertEqual(loaded, {"174": 0.10, "331": 0.90})
        for value in loaded.values():
            self.assertIsInstance(value, float)

    def test_fold_column_must_still_be_present_in_the_artifact(self) -> None:
        path = self._write(self.TABLE.replace("\touter_fold", "").replace("\t0\n", "\n").replace("\t3\n", "\n"))
        with self.assertRaises(PredictionTableError):
            load_predictions(
                path, arm="molecular", lineage_scope="all_lineage", endpoint="fibrosis"
            )

    def test_duplicate_out_of_fold_prediction_raises(self) -> None:
        path = self._write(
            self.TABLE + "174\tmolecular\tall_lineage\tfibrosis\t0.20\t1\n"
        )
        with self.assertRaises(PredictionTableError):
            load_predictions(
                path, arm="molecular", lineage_scope="all_lineage", endpoint="fibrosis"
            )


class UncertaintyTests(unittest.TestCase):
    def setUp(self) -> None:
        self.endpoint_rows = _rows(ENDPOINTS)
        self.fold_rows = _rows(FOLDS)
        self.fibrosis = endpoint_values("fibrosis", self.endpoint_rows)
        self.roster = registered_roster(
            "fibrosis", "pooled_all_donors", self.endpoint_rows, self.fold_rows
        )
        self.perfect = {donor: self.fibrosis[donor] for donor in self.roster.donor_ids}

    def test_bootstrap_resamples_donors_not_units(self) -> None:
        summary = pooled_donor_bootstrap(
            "auprc", self.perfect, self.fibrosis, self.roster, n_resamples=200
        )
        self.assertEqual(summary.n_donors, 37)
        self.assertEqual(summary.estimate, 1.0)
        self.assertGreaterEqual(summary.degenerate_resamples, 0)

    def test_degenerate_resamples_are_counted_not_dropped(self) -> None:
        summary = pooled_donor_bootstrap(
            "auprc", self.perfect, self.fibrosis, self.roster, n_resamples=200
        )
        self.assertIsInstance(summary.degenerate_resamples, int)
        self.assertIn("degenerate_resamples", summary.to_dict())

    def test_random_score_reference_is_above_prevalence(self) -> None:
        """AUPRC is upward biased at 15/37; prevalence is not the null."""

        reference = random_score_reference(
            "auprc", self.fibrosis, self.roster, n_draws=10_000
        )
        prevalence = sum(self.fibrosis.values()) / len(self.fibrosis)
        self.assertAlmostEqual(prevalence, 15 / 37, places=6)
        self.assertGreater(reference.mean, prevalence)
        self.assertAlmostEqual(reference.mean, 0.458, places=2)
        self.assertAlmostEqual(reference.percentile_95, 0.609, places=2)

    def test_steatosis_reference_matches_the_power_statement(self) -> None:
        steatosis = endpoint_values("steatosis", self.endpoint_rows)
        roster = registered_roster(
            "steatosis", "pooled_all_donors", self.endpoint_rows, self.fold_rows
        )
        reference = random_score_reference(
            "spearman", steatosis, roster, n_draws=10_000
        )
        self.assertAlmostEqual(reference.sd, 0.166, places=2)
        self.assertAlmostEqual(reference.percentile_95, 0.324, places=2)

    def test_permutation_null_inherits_the_prediction_tie_structure(self) -> None:
        """A tied scorer has a tighter null than a continuous one.

        ``self.perfect`` emits only two distinct values, so permuting it cannot
        reach the spread a continuous scorer reaches.  The two nulls answer
        different questions and neither substitutes for the other.
        """

        tied = permutation_null(
            "auprc", self.perfect, self.fibrosis, self.roster, n_permutations=10_000
        )
        continuous = random_score_reference(
            "auprc", self.fibrosis, self.roster, n_draws=10_000
        )
        self.assertLess(tied.percentile_95, continuous.percentile_95)
        self.assertLess(tied.sd, continuous.sd)


class FamilyTests(unittest.TestCase):
    FROZEN_PRIMARY = (
        "cholangiocyte",
        "fibroblast",
        "hepatocyte",
        "macrophage",
        "t_cell",
    )

    def test_family_is_twelve_molecular_tests(self) -> None:
        family = confirmatory_family(self.FROZEN_PRIMARY)
        self.assertEqual(len(family), 12)
        self.assertEqual({arm for arm, _, _ in family}, {"molecular"})

    def test_secondary_lineages_are_absent_from_the_family(self) -> None:
        family = confirmatory_family(self.FROZEN_PRIMARY)
        scopes = {scope for _, scope, _ in family}
        self.assertNotIn("endothelial_cell", scopes)
        self.assertNotIn("b_cell", scopes)

    def test_only_the_pooled_scope_is_confirmatory(self) -> None:
        self.assertEqual(scoring.CONFIRMATORY_SCOPES, frozenset({"pooled_all_donors"}))
        for scope in scoring.SCOPES - scoring.CONFIRMATORY_SCOPES:
            self.assertTrue(
                scope.startswith("adult_only") or scope.startswith("leave_one_well_out")
            )


if __name__ == "__main__":
    unittest.main()


class NullSidednessTests(unittest.TestCase):
    """Two nulls in one record must never share a name and differ in side.

    The frozen scores artifact carried
    ``continuous_random_score_reference.null_percentile_95`` = 0.3195, a
    two-sided |rho|, beside ``label_permutation_null.null_percentile_95`` =
    0.2711, a one-sided signed rho, under the same field name. Comparing an
    observed |rho| against the second runs a ten percent test wearing a five
    percent label. Average precision is one sided, so the fibrosis nulls agreed
    and hid it.
    """

    def setUp(self) -> None:
        self.endpoint_rows = _rows(ENDPOINTS)
        self.fold_rows = _rows(FOLDS)
        self.steatosis = endpoint_values("steatosis", self.endpoint_rows)
        self.fibrosis = endpoint_values("fibrosis", self.endpoint_rows)

    def _nulls(self, kind: str, endpoint: str, outcomes: dict[str, float]):
        roster = registered_roster(
            endpoint, "pooled_all_donors", self.endpoint_rows, self.fold_rows
        )
        predictions = {
            donor: float(index)
            for index, donor in enumerate(sorted(roster.donor_ids))
        }
        return (
            permutation_null(kind, predictions, outcomes, roster, n_permutations=6_000),
            random_score_reference(kind, outcomes, roster, n_draws=6_000),
        )

    def test_rank_nulls_share_one_side(self) -> None:
        permuted, continuous = self._nulls("spearman", "steatosis", self.steatosis)
        self.assertEqual(permuted.percentile_scale, "absolute_two_sided")
        self.assertEqual(continuous.percentile_scale, "absolute_two_sided")
        self.assertAlmostEqual(
            permuted.percentile_95, continuous.percentile_95, delta=0.02
        )

    def test_rank_percentile_matches_its_own_p_value_side(self) -> None:
        """A two-sided p-value beside a one-sided percentile is incoherent."""

        permuted, _ = self._nulls("spearman", "steatosis", self.steatosis)
        self.assertEqual(permuted.percentile_scale, "absolute_two_sided")
        self.assertNotAlmostEqual(permuted.percentile_95, 0.271, places=2)
        self.assertAlmostEqual(permuted.percentile_95, 0.32, delta=0.03)

    def test_average_precision_nulls_are_one_sided(self) -> None:
        permuted, continuous = self._nulls("auprc", "fibrosis", self.fibrosis)
        self.assertEqual(permuted.percentile_scale, "signed_one_sided")
        self.assertEqual(continuous.percentile_scale, "signed_one_sided")

    def test_every_emitted_null_declares_its_side(self) -> None:
        for kind, endpoint, outcomes in (
            ("spearman", "steatosis", self.steatosis),
            ("auprc", "fibrosis", self.fibrosis),
        ):
            for null in self._nulls(kind, endpoint, outcomes):
                payload = null.to_dict()
                self.assertIn("percentile_scale", payload)
                self.assertIn(
                    payload["percentile_scale"],
                    {"absolute_two_sided", "signed_one_sided"},
                )
