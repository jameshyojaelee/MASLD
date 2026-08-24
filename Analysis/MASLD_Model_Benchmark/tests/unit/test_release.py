from __future__ import annotations

from datetime import datetime
import json
import os
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

from masld_bench.artifacts import (
    canonical_hash,
    freeze_tree,
    sha256_file,
    write_json_exclusive,
    write_text_exclusive,
)
from masld_bench.firewall import authorize_one_time_outcome_join
from masld_bench.release import (
    SPEC_SCHEMA_VERSION,
    TIMESTAMP_SCHEMA_VERSION,
    ReleaseError,
    _TIMESTAMP_PROVIDER_VERIFIERS,
    _identity_sets,
    stage_release,
    verify_release,
)
from masld_bench.selection import verify_selection_lock
from masld_bench.tournament import verify_terminal_evaluation_authorization
from tests.unit._sealed_fixtures import (
    SELECTED_MODEL,
    TASK_IDS,
    commit_prediction_set,
    create_joint_outcome_bundle,
    create_scientific_selection,
    digest,
    freeze_power_set,
    freeze_terminal_authorization,
)


def _make_writable(root: Path) -> None:
    if not root.exists():
        return
    for path in (root, *root.rglob("*")):
        try:
            os.chmod(path, 0o755 if path.is_dir() else 0o644)
        except OSError:
            pass


def _fixture_timestamp_verifier(**kwargs: object) -> datetime:
    proof_path = Path(str(kwargs["proof_path"]))
    proof_sha256 = str(kwargs["proof_sha256"])
    if sha256_file(proof_path) != proof_sha256:
        raise ReleaseError("fixture provider proof hash mismatch")
    receipt_id = str(kwargs["provider_receipt_id"])
    bound_identity_sha256 = str(kwargs["bound_identity_sha256"])
    locator = str(kwargs["external_locator"])
    expected = f"{receipt_id}\n{bound_identity_sha256}\n{locator}\n"
    if proof_path.read_text(encoding="utf-8") != expected:
        raise ReleaseError("fixture provider proof payload mismatch")
    if receipt_id.startswith("provider-pre_unblind-"):
        return datetime.fromisoformat("2026-01-01T00:00:00+00:00")
    if receipt_id.startswith("provider-terminal-"):
        return datetime.fromisoformat("2026-01-02T00:00:00+00:00")
    raise ReleaseError("fixture provider receipt ID is unknown")


class OpenReleaseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.provider_patch = patch.dict(
            _TIMESTAMP_PROVIDER_VERIFIERS,
            {"open timestamp authority": _fixture_timestamp_verifier},
            clear=True,
        )
        cls.provider_patch.start()
        cls.temporary = tempfile.TemporaryDirectory()
        cls.root = Path(cls.temporary.name)
        scientific = create_scientific_selection(cls.root)
        cls.candidate = scientific.selection_candidate_ledger_dir
        cls.selection = scientific.selection_lock_dir
        cls.lock = verify_selection_lock(cls.selection)
        state = cls.root / "evaluator"
        cls.commits = commit_prediction_set(
            cls.root,
            selection_lock_dir=cls.selection,
            evaluator_state_dir=state,
        )
        cls.powers = freeze_power_set(evaluator_state_dir=state, commits=cls.commits)
        cls.outcome = create_joint_outcome_bundle(cls.root)
        consumption_ledger = cls.root / "independent-consumption-ledger"
        consumption_ledger.mkdir()
        consumption = authorize_one_time_outcome_join(
            evaluator_state_dir=state,
            power_decision_sha256s=(
                power.power_decision_sha256 for power in cls.powers.values()
            ),
            sealed_outcome_bundle_dir=cls.outcome,
            consumption_ledger_dir=consumption_ledger,
        )
        cls.authorization_dir = freeze_terminal_authorization(
            cls.root,
            selection_lock_dir=cls.selection,
            evaluator_state_dir=state,
            commits=cls.commits,
            powers=cls.powers,
            sealed_outcome_bundle_dir=cls.outcome,
            outcome_consumption_path=Path(consumption.marker_path),
        )
        cls.candidate_manifest_sha = sha256_file(cls.candidate / "ARTIFACTS.json")
        cls.selection_manifest_sha = sha256_file(cls.selection / "ARTIFACTS.json")
        cls.selection_lock_sha = sha256_file(cls.selection / "selection_lock.json")

    @classmethod
    def tearDownClass(cls) -> None:
        _make_writable(cls.root)
        cls.temporary.cleanup()
        cls.provider_patch.stop()

    def _authorization(
        self,
        *,
        label: str = "default",
    ) -> tuple[dict[str, object], dict[str, object]]:
        payload = json.loads(
            (self.authorization_dir / "terminal_evaluation_authorization.json").read_text(
                encoding="utf-8"
            )
        )
        return payload, {
            "path": self.authorization_dir.as_posix(),
            "manifest_sha256": sha256_file(
                self.authorization_dir / "ARTIFACTS.json"
            ),
            "document_sha256": sha256_file(
                self.authorization_dir
                / "terminal_evaluation_authorization.json"
            ),
        }

    def _forged_authorization(
        self, *, label: str, mutate
    ) -> tuple[dict[str, object], dict[str, object]]:
        payload, _ = self._authorization()
        forged = json.loads(json.dumps(payload))
        mutate(forged)
        identity = dict(forged)
        identity.pop("authorization_id", None)
        forged["authorization_id"] = canonical_hash(identity)
        root = self.root / f"forged-authorization-{label}"
        root.mkdir()
        document = write_json_exclusive(
            root / "terminal_evaluation_authorization.json", forged
        )
        manifest_sha = freeze_tree(
            root,
            {
                "artifact_class": "terminal_evaluation_authorization",
                "authorization_id": forged["authorization_id"],
            },
        )
        return forged, {
            "path": root.as_posix(),
            "manifest_sha256": manifest_sha,
            "document_sha256": document.sha256,
        }

    def _receipt(
        self,
        *,
        phase: str,
        identity: dict[str, str],
        timestamp: str,
        label: str = "default",
        provider: str = "Open Timestamp Authority",
        locator: str = "https://timestamp.example.org/receipts/immutable",
    ) -> tuple[dict[str, object], dict[str, object]]:
        root = self.root / f"receipt-{phase}-{label}"
        root.mkdir()
        provider_receipt_id = f"provider-{phase}-{label}"
        bound_identity_sha256 = canonical_hash(identity)
        proof = write_text_exclusive(
            root / "provider.proof",
            f"{provider_receipt_id}\n{bound_identity_sha256}\n{locator}\n",
        )
        normalized: dict[str, object] = {
            "schema_version": TIMESTAMP_SCHEMA_VERSION,
            "phase": phase,
            "provider": provider,
            "provider_receipt_id": provider_receipt_id,
            "external_locator": locator,
            "timestamp_utc": timestamp,
            "bound_identity": identity,
            "bound_identity_sha256": bound_identity_sha256,
            "proof": {
                "path": "provider.proof",
                "sha256": proof.sha256,
                "size_bytes": proof.size_bytes,
            },
        }
        payload = {"receipt_id": canonical_hash(normalized), **normalized}
        document = write_json_exclusive(root / "timestamp_receipt.json", payload)
        manifest_sha = freeze_tree(
            root,
            {
                "artifact_class": "external_timestamp_receipt",
                "receipt_id": payload["receipt_id"],
            },
        )
        return payload, {
            "path": root.as_posix(),
            "manifest_sha256": manifest_sha,
            "document_sha256": document.sha256,
        }

    def _open_artifacts(self, *, prediction_header: str | None = None) -> list[dict[str, object]]:
        root = self.root / "open_artifacts"
        root.mkdir(exist_ok=True)
        checkpoint_hashes = "\n".join(
            str(item["checkpoint_sha256"]) for item in self.lock.task_decisions
        )
        content = {
            "source.py": "SEED = 17\n",
            "environment.lock": "python==3.11.3\n",
            "model_card.md": f"# Model card\nModel: {SELECTED_MODEL}\nResearch use only.\n",
            "predictions.tsv": (
                prediction_header
                or "unit_hash\ttask_id\tmodel_id\tprediction\n"
            )
            + f"{digest('release-unit')}\t{TASK_IDS[0]}\t{SELECTED_MODEL}\t0.5\n",
            "benchmark.tsv": "task_id\tmetric\tvalue\tn_units\ncell_state_mapping\tmacro_f1\t0.72\t30\n",
            "fetch.md": (
                f"# Fetch instructions\n{SELECTED_MODEL}\n{checkpoint_hashes}\n"
                "Model weights are not redistributed. Fetch or rebuild them from the locked source.\n"
            ),
        }
        classes = {
            "source.py": "source",
            "environment.lock": "environment",
            "model_card.md": "model_card",
            "predictions.tsv": "frozen_predictions",
            "benchmark.tsv": "benchmark_table",
            "fetch.md": "fetch_instructions",
        }
        records: list[dict[str, object]] = []
        for filename, text in content.items():
            path = root / filename
            if path.exists():
                path.unlink()
            record = write_text_exclusive(path, text)
            records.append(
                {
                    "class": classes[filename],
                    "path": path.resolve().as_posix(),
                    "sha256": record.sha256,
                    "size_bytes": record.size_bytes,
                    "destination": filename,
                }
            )
        return records

    def _spec(
        self,
        *,
        authorization_payload: dict[str, object] | None = None,
        authorization_binding: dict[str, object] | None = None,
        artifacts: list[dict[str, object]] | None = None,
        pre_identity_override: dict[str, str] | None = None,
        pre_provider: str = "Open Timestamp Authority",
        pre_locator: str = "https://timestamp.example.org/receipts/pre",
        label: str = "default",
    ) -> dict[str, object]:
        if authorization_payload is None or authorization_binding is None:
            authorization_payload, authorization_binding = self._authorization(label=label)
        identity_sets = _identity_sets(authorization_payload)
        pre_identity = pre_identity_override or {
            "selection_lock_id": self.lock.lock_id,
            "selection_lock_manifest_sha256": self.selection_manifest_sha,
            "prediction_set_sha256": identity_sets["prediction_set_sha256"],
            "power_set_sha256": identity_sets["power_set_sha256"],
        }
        _, pre_binding = self._receipt(
            phase="pre_unblind",
            identity=pre_identity,
            timestamp="2026-01-01T00:00:00Z",
            label=label,
            provider=pre_provider,
            locator=pre_locator,
        )
        terminal_identity = {
            "authorization_id": str(authorization_payload["authorization_id"]),
            "authorization_manifest_sha256": str(authorization_binding["manifest_sha256"]),
            "outcome_set_sha256": identity_sets["outcome_set_sha256"],
            "terminal_decisions_sha256": identity_sets["terminal_decisions_sha256"],
        }
        _, terminal_binding = self._receipt(
            phase="terminal",
            identity=terminal_identity,
            timestamp="2026-01-02T00:00:00Z",
            label=label,
        )
        return {
            "schema_version": SPEC_SCHEMA_VERSION,
            "release_id": f"masld-open-{label}",
            "candidate": {
                "path": self.candidate.resolve().as_posix(),
                "manifest_sha256": self.candidate_manifest_sha,
            },
            "selection_lock": {
                "path": self.selection.resolve().as_posix(),
                "manifest_sha256": self.selection_manifest_sha,
                "lock_sha256": self.selection_lock_sha,
            },
            "terminal_authorization": authorization_binding,
            "receipts": {
                "pre_unblind": pre_binding,
                "terminal": terminal_binding,
            },
            "artifacts": artifacts or self._open_artifacts(),
        }

    def _write_spec(self, spec: dict[str, object], label: str) -> Path:
        path = self.root / f"spec-{label}.json"
        write_json_exclusive(path, spec)
        return path

    def _stage(self, spec: dict[str, object], label: str) -> Path:
        return stage_release(
            spec_path=self._write_spec(spec, label),
            output_root=self.root / "releases",
        )

    def test_stages_verified_terminal_failure_release(self) -> None:
        release = self._stage(self._spec(label="terminal-failure"), "terminal-failure")
        manifest = verify_release(release)
        self.assertEqual(
            manifest["publication_disposition"], "negative_result_no_champion"
        )
        self.assertEqual(manifest["champions"], [])
        self.assertFalse((release / "release_manifest.json").stat().st_mode & 0o222)

    def test_all_terminal_failures_produce_valid_zero_champion_release(self) -> None:
        authorization, binding = self._authorization(label="negative")
        spec = self._spec(
            authorization_payload=authorization,
            authorization_binding=binding,
            label="negative",
        )
        release = self._stage(spec, "negative")
        manifest = verify_release(release)
        self.assertEqual(manifest["champions"], [])
        self.assertEqual(
            manifest["publication_disposition"], "negative_result_no_champion"
        )

    def test_rejects_selection_lock_document_hash_mismatch(self) -> None:
        spec = self._spec(label="bad-lock")
        spec["selection_lock"]["lock_sha256"] = digest("wrong-lock")  # type: ignore[index]
        with self.assertRaisesRegex(ReleaseError, "SelectionLock document"):
            self._stage(spec, "bad-lock")

    def test_rejects_authorization_that_omits_locked_task(self) -> None:
        authorization, binding = self._forged_authorization(
            label="missing-task",
            mutate=lambda payload: payload["decisions"].pop(),
        )
        spec = self._spec(
            authorization_payload=authorization,
            authorization_binding=binding,
            label="missing-task",
        )
        with self.assertRaisesRegex(ReleaseError, "recursive verification"):
            self._stage(spec, "missing-task")

    def test_rejects_promoted_unsealed_decision(self) -> None:
        def mutate(payload: dict[str, object]) -> None:
            first = payload["decisions"][0]
            first["sealed"] = False
            first["promoted"] = True
            first["gate_passed"] = True

        authorization, binding = self._forged_authorization(
            label="unsealed", mutate=mutate
        )
        spec = self._spec(
            authorization_payload=authorization,
            authorization_binding=binding,
            label="unsealed",
        )
        with self.assertRaisesRegex(ReleaseError, "recursive verification"):
            self._stage(spec, "unsealed")

    def test_rejects_pre_unblind_receipt_bound_to_wrong_prediction_set(self) -> None:
        authorization, binding = self._authorization(label="wrong-pre")
        sets = _identity_sets(authorization)
        wrong = {
            "selection_lock_id": self.lock.lock_id,
            "selection_lock_manifest_sha256": self.selection_manifest_sha,
            "prediction_set_sha256": digest("wrong-predictions"),
            "power_set_sha256": sets["power_set_sha256"],
        }
        spec = self._spec(
            authorization_payload=authorization,
            authorization_binding=binding,
            pre_identity_override=wrong,
            label="wrong-pre",
        )
        with self.assertRaisesRegex(ReleaseError, "binds the wrong identities"):
            self._stage(spec, "wrong-pre")

    def test_rejects_local_timestamp_provider_and_locator(self) -> None:
        spec = self._spec(
            label="local-receipt",
            pre_provider="local",
            pre_locator="https://localhost/receipt",
        )
        with self.assertRaisesRegex(ReleaseError, "independent and non-local"):
            self._stage(spec, "local-receipt")

    def test_rejects_public_url_without_a_provider_specific_verifier(self) -> None:
        spec = self._spec(
            label="unknown-provider",
            pre_provider="Unregistered Timestamp Vendor",
            pre_locator="https://timestamp.example.org/receipts/pre",
        )
        with self.assertRaisesRegex(ReleaseError, "provider-specific.*fail-closed"):
            self._stage(spec, "unknown-provider")

    def test_rejects_checkpoint_or_weight_artifact_class(self) -> None:
        artifacts = self._open_artifacts()
        artifacts[0]["class"] = "weights"
        spec = self._spec(artifacts=artifacts, label="weights")
        with self.assertRaisesRegex(ReleaseError, "prohibited or unknown class"):
            self._stage(spec, "weights")

    def test_rejects_exact_protected_file_as_open_release_source(self) -> None:
        verified = verify_terminal_evaluation_authorization(
            self.authorization_dir, reverify_sources=True
        )
        protected = verified["_protected_source_roots"]
        protected_file = next(Path(path) for path in protected if Path(path).is_file())
        artifacts = self._open_artifacts()
        artifacts[0].update(
            {
                "path": protected_file.as_posix(),
                "sha256": sha256_file(protected_file),
                "size_bytes": protected_file.stat().st_size,
            }
        )
        spec = self._spec(artifacts=artifacts, label="protected-file")
        with self.assertRaisesRegex(ReleaseError, "protected frozen tree"):
            self._stage(spec, "protected-file")

    def test_rejects_unit_level_fields_in_frozen_predictions(self) -> None:
        artifacts = self._open_artifacts(
            prediction_header="unit_hash\tdonor_id\tmodel_id\tprediction\n"
        )
        spec = self._spec(artifacts=artifacts, label="donor-data")
        with self.assertRaisesRegex(ReleaseError, "forbidden/non-output fields"):
            self._stage(spec, "donor-data")

    def test_verify_rejects_writable_release_tree(self) -> None:
        release = self._stage(self._spec(label="writable"), "writable")
        artifact = next((release / "artifacts").rglob("source.py"))
        os.chmod(artifact, 0o644)
        with self.assertRaisesRegex(ReleaseError, "not read-only"):
            verify_release(release)


if __name__ == "__main__":
    unittest.main()
