"""Guard the three mpralegnet sites that are deliberately NOT changed.

MPRALegNet is terminally blocked from scoring.  That creates pressure to
harmonise every nearby field to agree with the block.  Three sites assert the
model is a roster "candidate", and all three are correct: the defect is
provenance equivalence, not terms, and no roster status value expresses
"terminally blocked on provenance".  Two of the three are working drift guards
keyed to the roster value.

These tests fail if someone edits those sites to agree with the terminal
disposition.  The reasoning is recorded in
config/artifacts/incidents/mpralegnet_deliberate_non_changes_20260825.json.
"""

from __future__ import annotations

import json
from pathlib import Path
import tomllib
import unittest


ROOT = Path(__file__).resolve().parents[2]
RECORD = (
    ROOT
    / "config/artifacts/incidents/mpralegnet_deliberate_non_changes_20260825.json"
)


class MpralegnetDeliberateNonChangesTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.record = json.loads(RECORD.read_text(encoding="utf-8"))

    def test_record_declares_exactly_the_three_non_changed_sites(self) -> None:
        sites = {entry["site"] for entry in self.record["non_changes"]}
        self.assertEqual(sites, {"D", "E", "F"})
        self.assertEqual(self.record["files_modified"], 0)
        self.assertEqual(self.record["pins_touched"], 0)
        self.assertIs(self.record["base_registry_unchanged"], True)

    def test_site_D_roster_status_is_still_candidate_and_admission_blocking(self) -> None:
        """The defect is provenance, not terms; blocked_terms would be false."""
        with (ROOT / "config/models/regulatory_local.toml").open("rb") as handle:
            family = tomllib.load(handle)
        entry = next(
            model
            for model in family["models"]
            if model["model_id"] == "mpralegnet"
        )
        self.assertEqual(entry["status"], "candidate")
        self.assertIs(entry["admission_blocking"], True)
        # The licence is what makes blocked_terms the wrong label.
        self.assertEqual(entry["license_status"], "code_MIT_weights_MIT")

    def test_site_E_contract_test_drift_guard_is_intact(self) -> None:
        """A guard on the roster identity, not a claim that the model is scorable."""
        source = (ROOT / "tests/contract/test_registry_files.py").read_text(
            encoding="utf-8"
        )
        self.assertIn('"mpralegnet": (', source)
        marker = source.index('"mpralegnet": (')
        window = source[marker : marker + 320]
        self.assertIn('"pretrained_adapter"', window)
        self.assertIn('"candidate"', window)
        self.assertIn('"target_label_unexposed"', window)

    def test_site_F_row_universe_drift_guard_is_intact(self) -> None:
        """Second independent guard, keyed to the same roster value."""
        source = (
            ROOT / "scripts/build_gse281364_five_seed_row_universe.py"
        ).read_text(encoding="utf-8")
        self.assertIn(
            '"mpralegnet": ("candidate", "target_label_unexposed"),',
            source,
        )

    def test_terminal_disposition_still_says_scoring_is_unauthorised(self) -> None:
        """The non-changes are only correct while the block itself stands."""
        disposition = json.loads(
            (
                ROOT
                / "config/artifacts/models/mpralegnet/terminal_disposition_20260825.json"
            ).read_text(encoding="utf-8")
        )
        self.assertIs(
            disposition["binding_consequences"]["scoring_authorized"], False
        )
        self.assertIs(
            disposition["terminal_receipt"]["provenance_equivalence_pass"], False
        )

    def test_record_admits_the_block_is_currently_unenforceable(self) -> None:
        """The block is on paper only until C1 lands; the record must say so.

        This is the guard against the correction reading as closure.  If the
        capability registry ever stops granting mpralegnet an endpoint, this
        test should be revisited together with the record -- not deleted to
        make the suite quiet.
        """
        block = self.record["known_unenforceable_block"]
        self.assertIn("SCORABLE IN PRACTICE", block["statement"])
        self.assertIn("C1", block["pending_work"])
        self.assertTrue(block["do_not_read_this_record_as_closure"])

        # The claim must stay true: the gating registry still grants endpoints.
        import tomllib

        with (
            ROOT / "config/evaluation/variant_to_regulation_capabilities.toml"
        ).open("rb") as handle:
            capabilities = tomllib.load(handle)
        entry = next(
            model
            for model in capabilities["models"]
            if model["model_id"] == "mpralegnet"
        )
        self.assertTrue(
            entry["allowed_endpoints"],
            "capability registry no longer grants endpoints; the "
            "unenforceable-block record is now stale and must be revisited",
        )


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
