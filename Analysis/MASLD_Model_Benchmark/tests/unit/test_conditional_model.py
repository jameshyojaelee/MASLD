from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

from masld_bench.artifacts import canonical_hash, freeze_tree, write_json_exclusive
from masld_bench.conditional_model import (
    STACKING_EVIDENCE_SCHEMA_VERSION,
    ComplementarityEvidence,
    ConditionalModelError,
    build_context_model_spec,
)
from masld_bench.hashing import sha256_file


class ConditionalModelTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.evidence_counter = 0

    def tearDown(self) -> None:
        self.temporary.cleanup()

    @staticmethod
    def candidate(
        model_id: str,
        *,
        development_score: float,
        compute_cost: float = 1.0,
        exposure_status: str = "clean_declared",
        within_one_se: bool = True,
    ) -> dict[str, object]:
        decision = {
            "schema_version": "masld-bench-open-champion-license-decision-v1",
            "eligible": True,
            "exposure_status": exposure_status,
        }
        return {
            "candidate_id": canonical_hash({"candidate": model_id}),
            "model_id": model_id,
            "adaptation_regime": "common_lane",
            "candidate_configuration_sha256": canonical_hash(
                {"configuration": model_id}
            ),
            "selection_artifact_ensemble_sha256s": {
                "checkpoint_sha256": canonical_hash(
                    {"checkpoint": model_id}
                ),
                "preprocessing_sha256": canonical_hash(
                    {"preprocessing": model_id}
                ),
                "task_head_sha256": canonical_hash({"head": model_id}),
                "calibration_sha256": canonical_hash(
                    {"calibration": model_id}
                ),
            },
            "admitted": True,
            "within_one_se": within_one_se,
            "development_score": development_score,
            "compute_cost": compute_cost,
            "open_champion_eligible": True,
            "open_champion_license_decision": decision,
            "open_champion_license_decision_sha256": canonical_hash(decision),
        }

    def evidence(self, *, residual_correlation: float = 0.70) -> ComplementarityEvidence:
        self.evidence_counter += 1
        base = self.root / f"evidence-{self.evidence_counter}"
        source_bindings = []
        for role, count in (
            ("development_gain_bundle", 2),
            ("development_residual_bundle", 2),
            ("development_shortlist", 2),
        ):
            artifact_class, filename = {
                "development_gain_bundle": (
                    "development_stack_gain_bundle",
                    "development_stack_gain_bundle.json",
                ),
                "development_residual_bundle": (
                    "development_residual_bundle",
                    "development_residual_bundle.json",
                ),
                "development_shortlist": (
                    "development_shortlist",
                    "development_shortlist.json",
                ),
            }[role]
            for index in range(count):
                source = base / f"source-{role}-{index}"
                source.mkdir(parents=True)
                write_json_exclusive(
                    source / filename,
                    {"fixture": True, "role": role, "index": index},
                )
                freeze_tree(source, {"artifact_class": artifact_class})
                source_bindings.append(
                    {
                        "role": role,
                        "artifact_class": artifact_class,
                        "path": source.resolve().as_posix(),
                        "manifest_sha256": sha256_file(source / "ARTIFACTS.json"),
                        "document_filename": filename,
                        "document_sha256": sha256_file(source / filename),
                    }
                )
        source_bindings.sort(key=lambda item: (item["role"], item["path"]))
        derived = {
            "residual_correlation": residual_correlation,
            "relative_deviance_gains": {"profile_deviance": 0.06},
            "absolute_correlation_or_f1_gains": {"variant_correlation": 0.03},
            "qualifying_seed_count": 4,
            "evaluated_seed_count": 5,
            "study_count": 2,
            "cross_fitted": True,
            "nonnegative_stack": True,
            "best_open_models_compared": True,
        }
        identity = {
            "schema_version": STACKING_EVIDENCE_SCHEMA_VERSION,
            "derived_evidence": derived,
            "source_bindings": source_bindings,
            "created_before_outcome_unblind": True,
            "sealed_results_used": False,
        }
        evidence_id = canonical_hash(identity)
        evidence_root = base / "stacking-evidence"
        evidence_root.mkdir()
        write_json_exclusive(
            evidence_root / "stacking_evidence.json",
            {"evidence_id": evidence_id, **identity},
        )
        freeze_tree(
            evidence_root,
            {"artifact_class": "stacking_evidence", "evidence_id": evidence_id},
        )
        return ComplementarityEvidence.from_stacking_evidence(evidence_root)

    def test_direct_construction_can_never_authorize(self) -> None:
        """Only a frozen stacking-evidence artifact may authorize the trigger."""

        with self.assertRaises(TypeError):
            ComplementarityEvidence(
                residual_correlation=0.70,
                relative_deviance_gains={"profile_deviance": 0.06},
                absolute_correlation_or_f1_gains={"variant_correlation": 0.03},
                qualifying_seed_count=4,
                evaluated_seed_count=5,
                study_count=2,
                cross_fitted=True,
                nonnegative_stack=True,
                best_open_models_compared=True,
            )

    def test_mapping_cannot_impersonate_a_frozen_stacking_tree(self) -> None:
        with self.assertRaises(ConditionalModelError):
            ComplementarityEvidence.from_stacking_evidence(  # type: ignore[arg-type]
                {"schema_version": STACKING_EVIDENCE_SCHEMA_VERSION}
            )

    def test_from_stacking_evidence_rejects_a_foreign_artifact(self) -> None:
        foreign = self.root / "foreign"
        foreign.mkdir()
        write_json_exclusive(foreign / "foreign.json", {"foreign": True})
        freeze_tree(foreign, {"artifact_class": "foreign"})
        with self.assertRaises(ConditionalModelError):
            ComplementarityEvidence.from_stacking_evidence(foreign)

    def test_post_load_tampering_revokes_authorization(self) -> None:
        evidence = self.evidence()
        document = Path(evidence.stacking_evidence_binding["path"]) / (
            "stacking_evidence.json"
        )
        payload = json.loads(document.read_text(encoding="utf-8"))
        payload["derived_evidence"]["study_count"] = 99
        document.chmod(0o640)
        document.write_text(json.dumps(payload), encoding="utf-8")
        with self.assertRaisesRegex(ConditionalModelError, "stacking-evidence"):
            evidence.require_trigger()

    def test_trigger_is_fail_closed(self) -> None:
        with self.assertRaises(ConditionalModelError):
            self.evidence(residual_correlation=0.80).require_trigger()

    def test_authorized_spec_ranks_on_development_score_not_compute(self) -> None:
        """0.79 vs 0.80 is not a tie, so compute is never consulted here.

        The former name claimed this exercised a compute tie-break; it did not.
        The real tie behaviour is covered by the two tests below.
        """

        evidence = self.evidence()
        spec = build_context_model_spec(
            evidence=evidence,
            sequence_candidates=[
                self.candidate("large", development_score=0.79, compute_cost=10),
                self.candidate("small", development_score=0.80, compute_cost=2),
            ],
            cell_candidates=[
                self.candidate("cell", development_score=0.8)
            ],
        )
        self.assertEqual(spec["sequence_backbone"], "small")
        self.assertEqual(
            spec["sequence_candidate_binding"]["candidate_id"],
            canonical_hash({"candidate": "small"}),
        )
        self.assertEqual(
            spec["cell_candidate_binding"]["candidate_id"],
            canonical_hash({"candidate": "cell"}),
        )
        self.assertEqual(spec["hotspot_117_usage"], "read_only_post_hoc_projection_only")

    def test_exact_score_tie_without_a_measured_compute_basis_fails_closed(
        self,
    ) -> None:
        with self.assertRaisesRegex(ConditionalModelError, "compute basis"):
            build_context_model_spec(
                evidence=self.evidence(),
                sequence_candidates=[
                    self.candidate("large", development_score=0.80, compute_cost=10),
                    self.candidate("small", development_score=0.80, compute_cost=2),
                ],
                cell_candidates=[self.candidate("cell", development_score=0.8)],
            )

    def test_exact_score_tie_uses_compute_only_under_one_hash_bound_basis(
        self,
    ) -> None:
        basis = canonical_hash({"compute_basis": "fixture-node-hours-v1"})
        large = self.candidate("large", development_score=0.80, compute_cost=10)
        small = self.candidate("small", development_score=0.80, compute_cost=2)
        for candidate in (large, small):
            candidate["compute_cost_status"] = "measured"
            candidate["compute_basis_sha256"] = basis
        spec = build_context_model_spec(
            evidence=self.evidence(),
            sequence_candidates=[large, small],
            cell_candidates=[self.candidate("cell", development_score=0.8)],
        )
        self.assertEqual(spec["sequence_backbone"], "small")

        # One candidate measured on a different basis is not comparable.
        small["compute_basis_sha256"] = canonical_hash({"compute_basis": "other"})
        with self.assertRaisesRegex(ConditionalModelError, "compute basis"):
            build_context_model_spec(
                evidence=self.evidence(),
                sequence_candidates=[large, small],
                cell_candidates=[self.candidate("cell", development_score=0.8)],
            )

    def test_candidate_license_decision_is_hash_bound_and_cell_must_be_clean(self) -> None:
        tampered = self.candidate("sequence", development_score=0.8)
        tampered["open_champion_license_decision"]["eligible"] = False
        with self.assertRaisesRegex(ConditionalModelError, "no open long-range"):
            build_context_model_spec(
                evidence=self.evidence(),
                sequence_candidates=[tampered],
                cell_candidates=[self.candidate("cell", development_score=0.8)],
            )

        with self.assertRaisesRegex(ConditionalModelError, "no clean open cell-state"):
            build_context_model_spec(
                evidence=self.evidence(),
                sequence_candidates=[self.candidate("sequence", development_score=0.8)],
                cell_candidates=[
                    self.candidate(
                        "cell",
                        development_score=0.8,
                        exposure_status="target_label_unexposed",
                    )
                ],
            )


if __name__ == "__main__":
    unittest.main()
