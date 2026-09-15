"""Lock the frozen decisions in the donor-by-lineage phenotype campaign spec."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

_spec = importlib.util.spec_from_file_location(
    "freeze_campaign_spec",
    ROOT / "scripts/freeze_gse296875_phenotype_campaign_spec.py",
)
assert _spec is not None and _spec.loader is not None
campaign = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(campaign)


class FrozenPartitionTests(unittest.TestCase):
    def test_primary_lineages_are_the_frozen_membership_five(self) -> None:
        self.assertEqual(
            campaign.PRIMARY_LINEAGES,
            ("cholangiocyte", "fibroblast", "hepatocyte", "macrophage", "t_cell"),
        )

    def test_secondary_lineages_are_reported_but_not_gate_members(self) -> None:
        self.assertEqual(
            campaign.SECONDARY_LINEAGES, ("endothelial_cell", "b_cell")
        )
        family = campaign.confirmatory_family(campaign.PRIMARY_LINEAGES)
        scopes = {scope for _, scope, _ in family}
        for lineage in campaign.SECONDARY_LINEAGES:
            self.assertNotIn(lineage, scopes)

    def test_collapse_map_covers_every_author_label(self) -> None:
        self.assertEqual(len(campaign.COLLAPSE_MAP), 7)
        self.assertEqual(
            set(campaign.COLLAPSE_MAP.values()),
            set(campaign.PRIMARY_LINEAGES) | set(campaign.SECONDARY_LINEAGES),
        )


class BoundFixtureTests(unittest.TestCase):
    def test_every_bound_fixture_still_verifies(self) -> None:
        verified = campaign.verify_bound_fixtures(ROOT)
        self.assertEqual(len(verified), 6)
        self.assertEqual(
            verified["donor_folds"]["artifacts_sha256"],
            "e6bc9161acaab508c298c49036af1ad10a720949178a14cba42282933195d416",
        )

    def test_a_changed_fixture_is_rejected(self) -> None:
        original = dict(campaign.BOUND_FIXTURES)
        campaign.BOUND_FIXTURES["donor_folds"] = (
            original["donor_folds"][0],
            "0" * 64,
        )
        try:
            with self.assertRaises(campaign.CampaignSpecError):
                campaign.verify_bound_fixtures(ROOT)
        finally:
            campaign.BOUND_FIXTURES.update(original)


class MinimumCellTests(unittest.TestCase):
    def setUp(self) -> None:
        membership = ROOT / campaign.BOUND_FIXTURES["fragment_membership"][0]
        self.counts, self.donors = campaign.lineage_census(membership)

    def test_threshold_is_twenty_nuclei(self) -> None:
        self.assertEqual(campaign.MINIMUM_CELLS_PER_UNIT, 20)

    def test_primary_units_observed_under_the_frozen_partition(self) -> None:
        rows = campaign.threshold_sensitivity(self.counts, self.donors)
        selected = next(row for row in rows if row["selected"])
        self.assertEqual(selected["minimum_nuclei_per_unit"], 20)
        self.assertEqual(selected["primary_units_total"], 195)
        self.assertEqual(selected["primary_units_observed"], 170)
        self.assertEqual(selected["primary_units_masked"], 25)

    def test_the_thinnest_lineages_absorb_the_threshold(self) -> None:
        rows = campaign.census_rows(self.counts, self.donors)
        by_lineage = {row["lineage_id"]: row for row in rows}
        self.assertEqual(by_lineage["cholangiocyte"]["units_below_threshold"], 11)
        self.assertEqual(by_lineage["t_cell"]["units_below_threshold"], 10)
        self.assertEqual(by_lineage["hepatocyte"]["units_below_threshold"], 0)

    def test_frozen_primary_coverage(self) -> None:
        rows = campaign.census_rows(self.counts, self.donors)
        primary = sum(r["total_nuclei"] for r in rows if r["analysis_role"] == "primary")
        secondary = sum(
            r["total_nuclei"] for r in rows if r["analysis_role"] == "secondary"
        )
        self.assertEqual(primary, 59_989)
        self.assertEqual(secondary, 8_409)
        self.assertEqual(primary + secondary, 68_398)


class PowerStatementTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.power = campaign.power_statement(ROOT)

    def test_position_is_stated_before_fitting(self) -> None:
        self.assertEqual(self.power["donors"], 39)
        self.assertEqual(self.power["fibrosis"]["positives"], 15)
        self.assertEqual(self.power["fibrosis"]["negatives"], 22)
        self.assertEqual(self.power["steatosis"]["observed_donors"], 38)

    def test_prevalence_is_rejected_as_the_auprc_null(self) -> None:
        fibrosis = self.power["fibrosis"]
        self.assertIs(fibrosis["prevalence_is_not_the_null"], True)
        reference = fibrosis["random_score_reference"]
        self.assertGreater(reference["null_mean"], fibrosis["prevalence"])
        self.assertAlmostEqual(reference["null_mean"], 0.458, places=2)
        self.assertAlmostEqual(reference["null_percentile_95"], 0.609, places=2)

    def test_steatosis_threshold_is_recorded(self) -> None:
        reference = self.power["steatosis"]["random_score_reference"]
        self.assertAlmostEqual(reference["null_percentile_95"], 0.324, places=2)

    def test_steatosis_leverage_is_recorded(self) -> None:
        steatosis = self.power["steatosis"]
        self.assertEqual(steatosis["donors_at_zero_percent"], 10)
        self.assertEqual(steatosis["donors_at_or_below_ten_percent"], 29)
        self.assertEqual(steatosis["donors_at_maximum"], 2)
        self.assertEqual(steatosis["maximum_percent"], 85.0)
        self.assertIn("leverage", steatosis["leverage_note"])

    def test_a_negative_is_an_acceptable_outcome(self) -> None:
        self.assertIn("negative", self.power["acceptable_outcome"])
        self.assertIn("publishable", self.power["acceptable_outcome"])

    def test_ambient_confounding_is_named_as_undetectable(self) -> None:
        self.assertTrue(
            any("ambient" in item for item in self.power["cannot_detect"])
        )


class SpecContentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.spec, cls.census, cls.sensitivity = campaign.build_spec(ROOT)

    def test_race_is_never_fitted(self) -> None:
        self.assertIs(self.spec["arms"]["race_is_not_fitted"], True)
        self.assertNotIn("reported_race", self.spec["arms"]["metadata"]["features"])
        self.assertEqual(
            self.spec["arms"]["metadata"]["features"],
            ["age_in_yr", "reported_sex", "BMI"],
        )

    def test_family_is_twelve_molecular_tests(self) -> None:
        self.assertEqual(self.spec["multiple_testing"]["family_size"], 12)
        arms = {entry["arm"] for entry in self.spec["multiple_testing"]["family"]}
        self.assertEqual(arms, {"molecular"})

    def test_bootstrap_contract_is_satisfied_not_deviated_from(self) -> None:
        uncertainty = self.spec["uncertainty"]
        self.assertIs(uncertainty["contract_satisfied"], True)
        self.assertEqual(uncertainty["contract_field"], "paired_donor_cluster_bootstrap")
        self.assertEqual(uncertainty["resampling_unit"], "donor")
        self.assertEqual(uncertainty["replicates"], 10_000)

    def test_campaign_cannot_confirm_a_champion(self) -> None:
        self.assertIs(self.spec["champion_eligible"], False)
        self.assertEqual(self.spec["role"], "secondary_development_only")
        self.assertIn("champion confirmation", self.spec["claims"]["forbidden"])

    def test_hotspot_programs_are_not_targets(self) -> None:
        hotspot = self.spec["hotspot_programs"]
        self.assertIs(hotspot["used_as_targets"], False)
        self.assertIs(hotspot["used_as_selected_features"], False)
        self.assertIs(hotspot["used_as_tuning_endpoints"], False)

    def test_sensitivity_table_is_auditable(self) -> None:
        thresholds = [row["minimum_nuclei_per_unit"] for row in self.sensitivity]
        self.assertEqual(thresholds, [10, 20, 25, 30, 50])
        self.assertEqual(sum(row["selected"] for row in self.sensitivity), 1)

    def test_corgi_artifact_is_bound_as_a_cross_check_not_a_source(self) -> None:
        self.assertIn("corgi_context_counts_cross_check", self.spec["bound_fixtures"])
        self.assertEqual(
            self.spec["fixture"]["kind"],
            "new_full_transcriptome_donor_by_lineage_pseudobulk",
        )

    def test_environment_note_records_the_numpy_collision(self) -> None:
        self.assertIn("PYTHONNOUSERSITE", self.spec["environment_note"])
        self.assertIn("_ARRAY_API", self.spec["environment_note"])


if __name__ == "__main__":
    unittest.main()
