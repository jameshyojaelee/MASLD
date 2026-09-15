"""The census-binding scope must not become a way to skip census binding.

`test_every_campaign_binds_the_frozen_model_census` used to index
`campaign['selection']` on every config/campaigns/*.toml, which errored on the
five dataset-lane records that legitimately have no model census to bind. The
fix scopes the check by the model-dispatch surface. These tests exist to prove
the scoping did not open a hole: a campaign that dispatches models still cannot
avoid binding the census.
"""

from __future__ import annotations

import unittest

from tests.contract.test_registry_files import (
    CAMPAIGN_MODEL_DISPATCH_KEYS,
    campaign_dispatch_keys,
    references_registry_model_id,
)


def model_campaign(**overrides: object) -> dict:
    campaign = {key: "x" for key in CAMPAIGN_MODEL_DISPATCH_KEYS}
    campaign["selection"] = {"census_model_ids_sha256": "deadbeef"}
    campaign.update(overrides)
    return campaign


class DispatchSurfaceTests(unittest.TestCase):
    def test_a_full_model_campaign_is_recognised(self) -> None:
        self.assertEqual(
            campaign_dispatch_keys(model_campaign()),
            set(CAMPAIGN_MODEL_DISPATCH_KEYS),
        )

    def test_a_dataset_lane_record_declares_nothing(self) -> None:
        record = {
            "schema_version": "masld-bench-campaign-v1",
            "campaign_id": "gse999999_reprocessing_v1",
            "participants": 12,
            "models": [{"model_kind": "training_class_prior"}],
        }
        self.assertEqual(campaign_dispatch_keys(record), set())

    def test_a_model_campaign_missing_selection_is_not_mistaken_for_inert(self) -> None:
        """The whole point of the fix: this must still fail, not be skipped."""

        campaign = model_campaign()
        del campaign["selection"]
        declared = campaign_dispatch_keys(campaign)
        self.assertNotEqual(declared, set(), "would be skipped as dataset-lane")
        self.assertNotEqual(
            declared,
            set(CAMPAIGN_MODEL_DISPATCH_KEYS),
            "must be flagged as a partial declaration",
        )
        self.assertIn("selection", CAMPAIGN_MODEL_DISPATCH_KEYS - declared)

    def test_dropping_the_whole_surface_is_caught_by_the_model_id_rule(self) -> None:
        """Stripping every dispatch key must not launder a model campaign."""

        campaign = {
            "schema_version": "masld-bench-campaign-v1",
            "models": [{"model_id": "scgpt"}],
        }
        self.assertEqual(campaign_dispatch_keys(campaign), set())
        self.assertTrue(references_registry_model_id(campaign))


class RegistryModelReferenceTests(unittest.TestCase):
    def test_inline_model_kind_specs_are_not_registry_references(self) -> None:
        record = {
            "models": [
                {"model_kind": "training_class_prior", "seed_policy": "seed_0"},
                {"model_kind": "gene_rank_nearest_centroid"},
            ]
        }
        self.assertFalse(references_registry_model_id(record))

    def test_a_nested_model_id_is_found(self) -> None:
        self.assertTrue(
            references_registry_model_id({"a": {"b": [{"model_ids": ["uce"]}]}})
        )

    def test_plain_records_have_no_registry_reference(self) -> None:
        self.assertFalse(references_registry_model_id({"participants": 26}))
        self.assertFalse(references_registry_model_id([1, "two", None]))


if __name__ == "__main__":
    unittest.main()
