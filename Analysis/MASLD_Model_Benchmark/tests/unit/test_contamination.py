from __future__ import annotations

import unittest

from masld_bench.contamination import (
    ArtifactKind,
    ArtifactSignatureRecord,
    CheckpointAuditRecord,
    CohortIdentityRecord,
    ContaminationAuditError,
    ExposureConfidence,
    SaltedDonorFingerprintRecord,
    ScBaseCountSnapshot,
    artifact_duplicate_components,
    assess_checkpoint_champion_eligibility,
    donor_duplicate_components,
    resolve_cohort_families,
)
from masld_bench.contracts import ExposureState


def _digest(character: str) -> str:
    return character * 64


def _cohort_record(record_id: str, **aliases: tuple[str, ...]) -> CohortIdentityRecord:
    fields: dict[str, object] = {
        "record_id": record_id,
        "accession_aliases": (),
        "bioproject_aliases": (),
        "sra_run_aliases": (),
        "doi_aliases": (),
        "donor_aliases": (),
        "cohort_release_aliases": (),
    }
    fields.update(aliases)
    return CohortIdentityRecord(**fields)


def _fingerprint(
    subject: str,
    kind: str,
    fingerprint: str,
    *,
    salt: str = "f",
) -> SaltedDonorFingerprintRecord:
    return SaltedDonorFingerprintRecord(
        subject_record_id=subject,
        fingerprint_kind=kind,
        fingerprint_sha256=_digest(fingerprint),
        salt_commitment_sha256=_digest(salt),
    )


def _checkpoint(**overrides: object) -> CheckpointAuditRecord:
    fields: dict[str, object] = {
        "schema_version": "masld-bench-checkpoint-audit-v1",
        "model_id": "model-a",
        "checkpoint_revision": "release-1",
        "checkpoint_sha256": _digest("a"),
        "release_date": "2026-08-01",
        "training_cutoff": "2025-12-31",
        "declared_corpora": ("Human Cell Atlas release 1",),
        "cellxgene_uuids": (),
        "accession_aliases": ("GSE100000",),
        "corpus_declaration_complete": True,
        "exposure_state": ExposureState.CLEAN_DECLARED,
        "exposure_confidence": ExposureConfidence.HIGH,
    }
    fields.update(overrides)
    return CheckpointAuditRecord(**fields)


def _snapshots() -> tuple[ScBaseCountSnapshot, ...]:
    return (
        ScBaseCountSnapshot(
            snapshot_id="scbasecount-2025-12",
            release_date="2025-12-01",
            corpus_aliases=("scBaseCount snapshot 2025-12",),
            accession_aliases=("GSE100000",),
            cellxgene_uuids=("11111111-1111-4111-8111-111111111111",),
        ),
        ScBaseCountSnapshot(
            snapshot_id="scbasecount-2026-02-target",
            release_date="2026-02-15",
            corpus_aliases=("scBaseCount snapshot 2026-02",),
            accession_aliases=("GSE289173", "GSE200000"),
            cellxgene_uuids=("22222222-2222-4222-8222-222222222222",),
        ),
        ScBaseCountSnapshot(
            snapshot_id="scbasecount-2026-06-later",
            release_date="2026-06-01",
            corpus_aliases=("scBaseCount snapshot 2026-06",),
            accession_aliases=("GSE300000",),
            cellxgene_uuids=("33333333-3333-4333-8333-333333333333",),
        ),
    )


class CohortFamilyAuditTests(unittest.TestCase):
    def test_alias_resolution_is_transitive_and_order_independent(self) -> None:
        records = [
            _cohort_record("r1", accession_aliases=("gse296875",)),
            _cohort_record(
                "r2",
                accession_aliases=("GSE296875",),
                bioproject_aliases=("prjna123456",),
            ),
            _cohort_record(
                "r3",
                bioproject_aliases=("PRJNA123456",),
                sra_run_aliases=("srr123",),
            ),
            _cohort_record(
                "r4",
                sra_run_aliases=("SRR123",),
                doi_aliases=("https://doi.org/10.1000/Atlas",),
            ),
            _cohort_record(
                "r5",
                doi_aliases=("doi:10.1000/atlas",),
                donor_aliases=("10.1000/atlas::donor-7",),
            ),
            _cohort_record(
                "r6",
                donor_aliases=("10.1000/ATLAS::DONOR-7",),
                cohort_release_aliases=("10.1000/atlas::release-v2",),
            ),
            _cohort_record(
                "r7",
                cohort_release_aliases=("10.1000/atlas::release-v2",),
            ),
            _cohort_record("other", accession_aliases=("GSE999999",)),
        ]
        forward = resolve_cohort_families(records)
        reverse = resolve_cohort_families(tuple(reversed(records)))
        self.assertEqual(forward, reverse)
        self.assertEqual(sorted(len(family.record_ids) for family in forward), [1, 7])
        connected = next(family for family in forward if len(family.record_ids) == 7)
        self.assertEqual(connected.record_ids, tuple(f"r{index}" for index in range(1, 8)))
        self.assertRegex(connected.cohort_family_id, r"^cohort-family-[0-9a-f]{64}$")

    def test_unqualified_donor_alias_and_duplicate_record_ids_fail_closed(self) -> None:
        with self.assertRaisesRegex(ContaminationAuditError, "namespace::value"):
            _cohort_record("r1", donor_aliases=("donor-1",))
        record = _cohort_record("r1", accession_aliases=("GSE1",))
        with self.assertRaisesRegex(ContaminationAuditError, "duplicate.*record_id"):
            resolve_cohort_families((record, record))


class DonorAndArtifactAuditTests(unittest.TestCase):
    def test_salted_fingerprints_form_deterministic_duplicate_components(self) -> None:
        records = (
            _fingerprint("source-a::d1", "expression", "1"),
            _fingerprint("source-a::d1", "genotype", "2"),
            _fingerprint("source-a::d1", "metadata", "3"),
            _fingerprint("source-a::d1", "sex", "5"),
            _fingerprint("source-b::d9", "expression", "1"),
            _fingerprint("source-b::d9", "genotype", "2"),
            _fingerprint("source-b::d9", "metadata", "3"),
            _fingerprint("source-b::d9", "sex", "5"),
            _fingerprint("source-c::d2", "expression", "4"),
        )
        forward = donor_duplicate_components(records)
        reverse = donor_duplicate_components(tuple(reversed(records)))
        self.assertEqual(forward, reverse)
        self.assertEqual(len(forward), 1)
        self.assertEqual(
            forward[0].subject_record_ids,
            ("source-a::d1", "source-b::d9"),
        )
        self.assertEqual(
            tuple(kind.value for kind in forward[0].matching_fingerprint_kinds),
            ("expression", "genotype", "metadata", "sex"),
        )

    def test_expression_and_metadata_may_differ_for_the_same_donor(self) -> None:
        records = (
            _fingerprint("source-a::liver", "genotype", "1"),
            _fingerprint("source-a::liver", "sex", "2"),
            _fingerprint("source-a::liver", "expression", "3"),
            _fingerprint("source-a::liver", "metadata", "4"),
            _fingerprint("source-b::blood", "genotype", "1"),
            _fingerprint("source-b::blood", "sex", "2"),
            _fingerprint("source-b::blood", "expression", "5"),
            _fingerprint("source-b::blood", "metadata", "6"),
        )
        components = donor_duplicate_components(records)
        self.assertEqual(len(components), 1)
        self.assertEqual(
            components[0].subject_record_ids,
            ("source-a::liver", "source-b::blood"),
        )
        self.assertEqual(
            tuple(kind.value for kind in components[0].matching_fingerprint_kinds),
            ("genotype", "sex"),
        )

    def test_sex_or_metadata_match_alone_never_declares_duplicate_donors(self) -> None:
        records = (
            _fingerprint("source-a::d1", "sex", "1"),
            _fingerprint("source-b::d9", "sex", "1"),
            _fingerprint("source-a::d1", "metadata", "2"),
            _fingerprint("source-b::d9", "metadata", "2"),
        )
        self.assertEqual(donor_duplicate_components(records), ())

    def test_raw_or_conflicting_fingerprints_are_rejected(self) -> None:
        payload = {
            "subject_record_id": "source-a::d1",
            "fingerprint_kind": "expression",
            "fingerprint_sha256": _digest("1"),
            "salt_commitment_sha256": _digest("f"),
            "scheme": "sha256_salted_v1",
            "raw_expression": [1.0, 2.0],
        }
        with self.assertRaisesRegex(ContaminationAuditError, "unknown field"):
            SaltedDonorFingerprintRecord.from_dict(payload)
        with self.assertRaisesRegex(ContaminationAuditError, "64-character"):
            _fingerprint("source-a::d1", "expression", "not-a-hex-character")
        with self.assertRaisesRegex(ContaminationAuditError, "same salt"):
            donor_duplicate_components(
                (
                    _fingerprint("source-a::d1", "expression", "1", salt="e"),
                    _fingerprint("source-b::d9", "expression", "1", salt="f"),
                )
            )
        conflicting = (
            _fingerprint("source-a::d1", "expression", "1"),
            _fingerprint("source-b::d9", "expression", "1"),
            _fingerprint("source-a::d1", "genotype", "2"),
            _fingerprint("source-b::d9", "genotype", "3"),
        )
        with self.assertRaisesRegex(ContaminationAuditError, "conflicting.*genotype"):
            donor_duplicate_components(conflicting)
        conflicting_sex = (
            _fingerprint("source-a::d1", "expression", "1"),
            _fingerprint("source-b::d9", "expression", "1"),
            _fingerprint("source-a::d1", "sex", "2"),
            _fingerprint("source-b::d9", "sex", "3"),
        )
        with self.assertRaisesRegex(ContaminationAuditError, "conflicting.*sex"):
            donor_duplicate_components(conflicting_sex)

    def test_exact_and_supplied_near_artifact_signatures_are_grouped(self) -> None:
        records = (
            ArtifactSignatureRecord("m1", ArtifactKind.MATRIX, _digest("1"), None),
            ArtifactSignatureRecord(
                "m2", ArtifactKind.MATRIX, _digest("1"), _digest("a")
            ),
            ArtifactSignatureRecord(
                "m3", ArtifactKind.MATRIX, _digest("2"), _digest("a")
            ),
            ArtifactSignatureRecord(
                "i1", ArtifactKind.IMAGE, _digest("3"), _digest("b")
            ),
            ArtifactSignatureRecord(
                "i2", ArtifactKind.IMAGE, _digest("4"), _digest("b")
            ),
        )
        forward = artifact_duplicate_components(records)
        reverse = artifact_duplicate_components(tuple(reversed(records)))
        self.assertEqual(forward, reverse)
        by_kind = {component.artifact_kind: component for component in forward}
        self.assertEqual(by_kind[ArtifactKind.MATRIX].artifact_ids, ("m1", "m2", "m3"))
        self.assertEqual(
            by_kind[ArtifactKind.MATRIX].evidence,
            ("exact_sha256", "supplied_near_signature_sha256"),
        )
        self.assertEqual(by_kind[ArtifactKind.IMAGE].artifact_ids, ("i1", "i2"))
        self.assertEqual(
            by_kind[ArtifactKind.IMAGE].evidence,
            ("supplied_near_signature_sha256",),
        )

    def test_artifact_type_conflicts_are_rejected(self) -> None:
        with self.assertRaisesRegex(ContaminationAuditError, "both matrix and image"):
            artifact_duplicate_components(
                (
                    ArtifactSignatureRecord(
                        "matrix", ArtifactKind.MATRIX, _digest("1"), None
                    ),
                    ArtifactSignatureRecord(
                        "image", ArtifactKind.IMAGE, _digest("1"), None
                    ),
                )
            )


class CheckpointExposureAuditTests(unittest.TestCase):
    def test_checkpoint_record_is_strict_exact_and_canonical(self) -> None:
        record = _checkpoint(
            declared_corpora=("Z corpus", "a CORPUS"),
            accession_aliases=("gse100001", "GSE100000"),
        )
        reordered = _checkpoint(
            declared_corpora=("A corpus", "z corpus"),
            accession_aliases=("GSE100000", "GSE100001"),
        )
        self.assertEqual(record.audit_id, reordered.audit_id)
        self.assertRegex(record.audit_id, r"^[0-9a-f]{64}$")
        with self.assertRaisesRegex(ContaminationAuditError, "unknown field"):
            CheckpointAuditRecord.from_dict(
                {
                    **{
                        field: getattr(record, field)
                        for field in CheckpointAuditRecord._FIELDS
                    },
                    "weights_url": "https://example.invalid/model.bin",
                }
            )
        with self.assertRaisesRegex(ContaminationAuditError, "must be exact"):
            _checkpoint(checkpoint_revision="UNRESOLVED")
        with self.assertRaisesRegex(ContaminationAuditError, "must not follow"):
            _checkpoint(training_cutoff="2026-09-01")

    def test_only_two_high_confidence_exposure_states_are_champion_eligible(self) -> None:
        snapshots = _snapshots()
        for state in (ExposureState.CLEAN_DECLARED, ExposureState.TARGET_LABEL_UNEXPOSED):
            with self.subTest(state=state.value):
                decision = assess_checkpoint_champion_eligibility(
                    _checkpoint(exposure_state=state), scbasecount_snapshots=snapshots
                )
                self.assertTrue(decision.eligible)
                self.assertEqual(decision.reasons, ())
        for state in (
            ExposureState.ENCODER_SEEN,
            ExposureState.CONTINUAL_SEEN,
            ExposureState.REFERENCE_ONLY,
            ExposureState.DOWNSTREAM_DEMO,
            ExposureState.UNKNOWN,
        ):
            with self.subTest(state=state.value):
                decision = assess_checkpoint_champion_eligibility(
                    _checkpoint(exposure_state=state), scbasecount_snapshots=snapshots
                )
                self.assertFalse(decision.eligible)
                self.assertIn(
                    f"exposure_state_not_champion_eligible:{state.value}",
                    decision.reasons,
                )

    def test_target_bearing_and_later_scbasecount_snapshots_are_ineligible(self) -> None:
        snapshots = _snapshots()
        for corpus, expected_snapshot in (
            ("scBaseCount snapshot 2026-02", "scbasecount-2026-02-target"),
            ("scBaseCount snapshot 2026-06", "scbasecount-2026-06-later"),
        ):
            with self.subTest(corpus=corpus):
                decision = assess_checkpoint_champion_eligibility(
                    _checkpoint(
                        declared_corpora=(corpus,),
                        training_cutoff="2026-06-30",
                    ),
                    scbasecount_snapshots=snapshots,
                )
                self.assertFalse(decision.eligible)
                self.assertIn(expected_snapshot, decision.forbidden_snapshot_ids)
                self.assertIn(
                    f"forbidden_scbasecount_snapshot:{expected_snapshot}",
                    decision.reasons,
                )

        uuid_decision = assess_checkpoint_champion_eligibility(
            _checkpoint(
                declared_corpora=("CELLxGENE census exact datasets",),
                cellxgene_uuids=("33333333-3333-4333-8333-333333333333",),
                training_cutoff="2026-06-30",
            ),
            scbasecount_snapshots=snapshots,
        )
        self.assertFalse(uuid_decision.eligible)
        self.assertIn("scbasecount-2026-06-later", uuid_decision.forbidden_snapshot_ids)

    def test_direct_target_unknown_corpus_and_incomplete_audits_fail_closed(self) -> None:
        snapshots = _snapshots()
        cases = (
            (
                _checkpoint(accession_aliases=("GSE289173",)),
                "target_accession_declared:GSE289173",
            ),
            (
                _checkpoint(declared_corpora=("unknown",)),
                "unknown_or_undeclared_corpus:unknown",
            ),
            (
                _checkpoint(corpus_declaration_complete=False),
                "corpus_declaration_incomplete",
            ),
            (
                _checkpoint(exposure_confidence=ExposureConfidence.MEDIUM),
                "exposure_confidence_not_high:medium",
            ),
            (
                _checkpoint(declared_corpora=("scBaseCount snapshot 2026-09",)),
                "unregistered_scbasecount_snapshot:scbasecount snapshot 2026-09",
            ),
            (
                _checkpoint(declared_corpora=("CELLxGENE Census",)),
                "cellxgene_corpus_without_exact_uuids",
            ),
        )
        for record, reason in cases:
            with self.subTest(reason=reason):
                decision = assess_checkpoint_champion_eligibility(
                    record, scbasecount_snapshots=snapshots
                )
                self.assertFalse(decision.eligible)
                self.assertIn(reason, decision.reasons)

        no_registry = assess_checkpoint_champion_eligibility(
            _checkpoint(), scbasecount_snapshots=()
        )
        self.assertFalse(no_registry.eligible)
        self.assertIn(
            "target_2026_scbasecount_snapshot_registry_missing", no_registry.reasons
        )

    def test_checkpoint_decision_is_independent_of_snapshot_input_order(self) -> None:
        snapshots = _snapshots()
        forward = assess_checkpoint_champion_eligibility(
            _checkpoint(), scbasecount_snapshots=snapshots
        )
        reverse = assess_checkpoint_champion_eligibility(
            _checkpoint(), scbasecount_snapshots=tuple(reversed(snapshots))
        )
        self.assertEqual(forward, reverse)


if __name__ == "__main__":
    unittest.main()
