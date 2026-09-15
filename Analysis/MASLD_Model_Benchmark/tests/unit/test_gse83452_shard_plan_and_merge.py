from __future__ import annotations

import csv
import json
from pathlib import Path
import tempfile
import unittest
import unittest.mock

from scripts.merge_gse83452_gpl16686_shards import (
    EXPECTED_PACKAGE_VERSIONS,
    ShardMergeError,
    digest_column,
    merge,
    read_matrix,
    verify_shard,
)
from scripts.plan_gse83452_summarization_shards import (
    ShardPlanError,
    plan_shards,
)

ACCESSIONS = [f"GSM{2203254 + i}" for i in range(231)]


class ShardPlanTests(unittest.TestCase):
    def test_plan_covers_every_array_exactly_once_as_owner(self) -> None:
        plan, sentinels = plan_shards(ACCESSIONS, shards=8, sentinels=3)
        owned = [a for bucket in plan for a in bucket if a not in sentinels]
        self.assertEqual(sorted(owned + list(sentinels)), sorted(ACCESSIONS))
        self.assertEqual(len(set(owned)), len(owned))

    def test_sentinels_appear_in_every_shard(self) -> None:
        plan, sentinels = plan_shards(ACCESSIONS, shards=8, sentinels=3)
        self.assertEqual(len(sentinels), 3)
        for bucket in plan:
            for sentinel in sentinels:
                self.assertIn(sentinel, bucket)

    def test_sentinels_are_spread_not_clustered(self) -> None:
        _, sentinels = plan_shards(ACCESSIONS, shards=8, sentinels=3)
        positions = sorted(ACCESSIONS.index(s) for s in sentinels)
        self.assertGreater(positions[-1] - positions[0], len(ACCESSIONS) // 2)

    def test_assignment_is_deterministic(self) -> None:
        first, _ = plan_shards(ACCESSIONS, shards=8, sentinels=3)
        second, _ = plan_shards(ACCESSIONS, shards=8, sentinels=3)
        self.assertEqual(first, second)

    def test_bad_shard_or_sentinel_counts_are_rejected(self) -> None:
        with self.assertRaisesRegex(ShardPlanError, "shard count"):
            plan_shards(ACCESSIONS, shards=0, sentinels=3)
        with self.assertRaisesRegex(ShardPlanError, "sentinel count"):
            plan_shards(ACCESSIONS, shards=8, sentinels=0)


def shard_receipt(shard: int, per_array: list[dict[str, object]], **overrides: object) -> dict[str, object]:
    payload: dict[str, object] = {
        "status": "pass_label_blind_single_array_shard",
        "series": "GSE83452",
        "platform_id": "GPL16686",
        "shard": f"{shard:02d}",
        "arrays_read_per_summarization_call": 1,
        "matrices_round_trip_bitwise": True,
        "label_blind": True,
        "samples_excluded": 0,
        "package_versions": dict(EXPECTED_PACKAGE_VERSIONS),
        "per_array": per_array,
        "all_sample_RMA_run": False,
        "across_array_quantile_normalization_run": False,
        "GEO_series_matrix_read": False,
        "labels_read": False,
        "model_training_activated": False,
        "sealed_outcomes_read": False,
    }
    payload.update(overrides)
    return payload


class VerifyShardTests(unittest.TestCase):
    def test_clean_shard_passes(self) -> None:
        verify_shard(shard_receipt(0, []))

    def test_batched_shard_is_rejected(self) -> None:
        with self.assertRaisesRegex(ShardMergeError, "more than one array per call"):
            verify_shard(shard_receipt(0, [], arrays_read_per_summarization_call=4))

    def test_shard_without_round_trip_assertion_is_rejected(self) -> None:
        with self.assertRaisesRegex(ShardMergeError, "bitwise round trip"):
            verify_shard(shard_receipt(0, [], matrices_round_trip_bitwise=False))

    def test_shard_that_read_labels_is_rejected(self) -> None:
        with self.assertRaisesRegex(ShardMergeError, "labels_read"):
            verify_shard(shard_receipt(0, [], labels_read=True))


class MergeTests(unittest.TestCase):
    """Two shards, one sentinel, three features, exercised end to end."""

    FEATURES = ["16657436", "16657440", "16657450"]
    GENES = ["ENSG00000000003", "ENSG00000000005"]

    def _write_shard(
        self, base: Path, shard: int, columns: dict[str, list[float]],
        genes: dict[str, list[float]], *, digest_override: dict[str, str] | None = None,
    ) -> Path:
        directory = base / f"shard_{shard:02d}"
        directory.mkdir()
        accessions = sorted(columns)
        for name, idcol, axis, source in (
            (f"shard_{shard:02d}_platform_feature_matrix.tsv", "platform_feature_id", self.FEATURES, columns),
            (f"shard_{shard:02d}_gene_matrix.tsv", "ensembl_gene_id", self.GENES, genes),
        ):
            with (directory / name).open("w", encoding="utf-8", newline="") as handle:
                handle.write("\t".join([idcol, *accessions]) + "\n")
                for index, feature in enumerate(axis):
                    handle.write("\t".join([feature, *["%.17g" % source[a][index] for a in accessions]]) + "\n")
        per_array = []
        for accession in accessions:
            digest = (digest_override or {}).get(accession) or digest_column(columns[accession])
            per_array.append({"sample_accession": accession, "array_digest": digest})
        (directory / f"shard_{shard:02d}_receipt.json").write_text(
            json.dumps(shard_receipt(shard, per_array)), encoding="utf-8"
        )
        return directory

    def _fixture(self, base: Path, *, sentinel_drift: bool = False, **kwargs: object) -> dict[str, object]:
        sentinel = "GSM0000000"
        s_col = [1.5, 2.25, 3.125]
        s_gene = [1.5, 2.25]
        a_cols = {"GSM0000001": [4.5, 5.25, 6.125], sentinel: list(s_col)}
        a_genes = {"GSM0000001": [4.5, 5.25], sentinel: list(s_gene)}
        drifted = [v + 1e-13 for v in s_col] if sentinel_drift else list(s_col)
        b_cols = {"GSM0000002": [7.5, 8.25, 9.125], sentinel: drifted}
        b_genes = {"GSM0000002": [7.5, 8.25], sentinel: list(s_gene)}
        d0 = self._write_shard(base, 0, a_cols, a_genes, **kwargs)
        d1 = self._write_shard(base, 1, b_cols, b_genes, **kwargs)
        return {"shard_dirs": [d0, d1], "plan": {"sentinels": [sentinel]}}

    def test_merge_produces_all_arrays_once(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            base = Path(value)
            args = self._fixture(base)
            with unittest.mock.patch(
                "scripts.merge_gse83452_gpl16686_shards.EXPECTED_ARRAYS", 3
            ), unittest.mock.patch(
                "scripts.merge_gse83452_gpl16686_shards.EXPECTED_PLATFORM_FEATURES", 3
            ):
                receipt = merge(output=base / "merged", **args)
            self.assertEqual(receipt["arrays_summarized"], 3)
            self.assertIs(receipt["sentinel_digests_agree_across_shards"], True)
            self.assertIs(receipt["cross_job_single_array_reproducibility_verified"], True)
            rows, cols, _ = read_matrix(
                base / "merged/gse83452_gpl16686_platform_feature_matrix.tsv", "platform_feature_id"
            )
            self.assertEqual(cols, ["GSM0000000", "GSM0000001", "GSM0000002"])
            self.assertEqual(rows, self.FEATURES)

    def test_sentinel_drift_between_shards_is_rejected(self) -> None:
        """If independent jobs disagree, the single-array guarantee is void."""
        with tempfile.TemporaryDirectory() as value:
            base = Path(value)
            args = self._fixture(base, sentinel_drift=True)
            with unittest.mock.patch(
                "scripts.merge_gse83452_gpl16686_shards.EXPECTED_ARRAYS", 3
            ), unittest.mock.patch(
                "scripts.merge_gse83452_gpl16686_shards.EXPECTED_PLATFORM_FEATURES", 3
            ):
                with self.assertRaisesRegex(ShardMergeError, "does not hold across jobs"):
                    merge(output=base / "merged", **args)

    def test_column_digest_mismatch_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            base = Path(value)
            args = self._fixture(base, digest_override={"GSM0000001": "f" * 64})
            with unittest.mock.patch(
                "scripts.merge_gse83452_gpl16686_shards.EXPECTED_ARRAYS", 3
            ), unittest.mock.patch(
                "scripts.merge_gse83452_gpl16686_shards.EXPECTED_PLATFORM_FEATURES", 3
            ):
                with self.assertRaisesRegex(ShardMergeError, "digest does not match"):
                    merge(output=base / "merged", **args)

    def test_wrong_total_array_count_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            base = Path(value)
            args = self._fixture(base)
            with unittest.mock.patch(
                "scripts.merge_gse83452_gpl16686_shards.EXPECTED_PLATFORM_FEATURES", 3
            ):
                with self.assertRaisesRegex(ShardMergeError, "not 231"):
                    merge(output=base / "merged", **args)

    def test_plan_without_sentinels_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            base = Path(value)
            args = self._fixture(base)
            args["plan"] = {"sentinels": []}
            with self.assertRaisesRegex(ShardMergeError, "no sentinel"):
                merge(output=base / "merged", **args)


class DigestTests(unittest.TestCase):
    def test_digest_matches_little_endian_doubles(self) -> None:
        import struct
        from hashlib import sha256 as _sha
        values = [1.5, -2.25, 3.125]
        self.assertEqual(
            digest_column(values), _sha(struct.pack("<3d", *values)).hexdigest()
        )

    def test_digest_is_sensitive_to_a_one_ulp_change(self) -> None:
        import math
        base = [1.5, 2.25, 3.125]
        moved = [math.nextafter(base[0], math.inf), 2.25, 3.125]
        self.assertNotEqual(digest_column(base), digest_column(moved))


if __name__ == "__main__":
    unittest.main()
