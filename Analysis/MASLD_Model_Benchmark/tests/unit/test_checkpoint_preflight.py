from __future__ import annotations

import hashlib
import json
from pathlib import Path
import pickle
import tempfile
import unittest

from masld_bench.artifacts import verify_frozen_tree
from masld_bench.checkpoint_preflight import (
    CheckpointPreflightError,
    load_restricted_plain_mapping,
    stage_geneformer_bundle,
)


def _artifact(path: Path, relative: str) -> dict[str, object]:
    payload = path.read_bytes()
    return {
        "path": relative,
        "sha256": hashlib.sha256(payload).hexdigest(),
        "size_bytes": len(payload),
    }


class RestrictedPickleTests(unittest.TestCase):
    def test_plain_dictionary_round_trips(self) -> None:
        value = {"ENSG000001": 2.5, "<pad>": 0, "alias": "ENSG000002"}
        self.assertEqual(
            load_restricted_plain_mapping(pickle.dumps(value, protocol=4)), value
        )

    def test_global_construction_and_nonfinite_values_are_rejected(self) -> None:
        with self.assertRaisesRegex(
            CheckpointPreflightError, "global/class construction"
        ):
            load_restricted_plain_mapping(pickle.dumps(eval, protocol=4))
        with self.assertRaisesRegex(CheckpointPreflightError, "non-finite"):
            load_restricted_plain_mapping(
                pickle.dumps({"ENSG000001": float("nan")}, protocol=4)
            )
        with self.assertRaisesRegex(CheckpointPreflightError, "trailing content"):
            load_restricted_plain_mapping(
                pickle.dumps({"ENSG000001": 1}, protocol=4) + b"hidden"
            )


class NumericScalarAllowlistTests(unittest.TestCase):
    """The pickle allowlist admits numeric scalars and nothing else.

    Upstream Geneformer stores gene medians and token IDs as numpy scalars, so
    the stager must decode them.  Widening find_class is a security boundary
    change, so these tests pin exactly what stays refused.
    """

    def test_numeric_numpy_scalars_are_decoded_to_plain_python(self) -> None:
        numpy = _require_numpy(self)
        payload = pickle.dumps(
            {
                "median": numpy.float64(2.5),
                "token": numpy.int16(7),
                "plain": 3,
            }
        )
        decoded = load_restricted_plain_mapping(payload)
        self.assertEqual(decoded, {"median": 2.5, "token": 7, "plain": 3})
        for value in decoded.values():
            self.assertIn(type(value), (int, float))
            self.assertNotIn("numpy", type(value).__module__)

    def test_an_object_dtype_array_is_still_refused(self) -> None:
        numpy = _require_numpy(self)

        class Detonator:
            def __reduce__(self):
                return (print, ("this must never run",))

        payload = pickle.dumps({"k": numpy.array([Detonator()], dtype=object)})
        with self.assertRaises(CheckpointPreflightError):
            load_restricted_plain_mapping(payload)

    def test_arbitrary_callables_and_reduce_objects_are_still_refused(self) -> None:
        import os

        class Detonator:
            def __reduce__(self):
                return (print, ("this must never run",))

        for payload in (
            pickle.dumps({"k": Detonator()}),
            pickle.dumps({"k": os.system}),
            pickle.dumps({"k": eval}),
        ):
            with self.assertRaises(CheckpointPreflightError):
                load_restricted_plain_mapping(payload)

    def test_non_numeric_numpy_dtypes_are_refused(self) -> None:
        numpy = _require_numpy(self)
        for value in (
            numpy.str_("text"),
            numpy.datetime64("2026-08-23"),
            numpy.complex128(1 + 2j),
        ):
            with self.assertRaises(CheckpointPreflightError):
                load_restricted_plain_mapping(pickle.dumps({"k": value}))

    def test_decoding_needs_no_numpy_import(self) -> None:
        """The security-critical path must not execute numpy code at all.

        The stager runs under the bare project interpreter, which has no numpy,
        and decoding an untrusted pickle should never import it.
        """

        import masld_bench.checkpoint_preflight as preflight

        source = Path(preflight.__file__).read_text(encoding="utf-8")
        decoder = source[source.index("_NUMERIC_DTYPE_FORMATS = {"):source.index(
            "class _PlainDataUnpickler"
        )]
        self.assertNotIn("import numpy", decoder)

    def test_struct_decode_matches_numpy_exactly(self) -> None:
        numpy = _require_numpy(self)
        for value in (
            numpy.float64(2.001186019549122),
            numpy.float32(1.5),
            numpy.int16(-7),
            numpy.int64(1 << 40),
            numpy.uint8(255),
        ):
            decoded = load_restricted_plain_mapping(pickle.dumps({"k": value}))["k"]
            self.assertEqual(decoded, value.item())
            self.assertIn(type(decoded), (int, float, bool))

    def test_an_empty_string_key_is_carried_verbatim(self) -> None:
        """Upstream V1 gene_name_id_dict_gc30M has '' -> ENSG00000285325."""

        decoded = load_restricted_plain_mapping(
            pickle.dumps({"": "ENSG00000285325", "MT-TF": "ENSG00000210049"})
        )
        self.assertEqual(decoded[""], "ENSG00000285325")

    def test_non_string_keys_are_still_refused(self) -> None:
        with self.assertRaises(CheckpointPreflightError):
            load_restricted_plain_mapping(pickle.dumps({1: "a"}))


def _require_numpy(case: unittest.TestCase):
    try:
        import numpy
    except ImportError:  # pragma: no cover - environment dependent
        case.skipTest("numpy is unavailable in this environment")
    return numpy


class GeneformerStagingTests(unittest.TestCase):
    def _fixture(self, root: Path) -> tuple[Path, Path, Path]:
        source = root / "source"
        output = root / "output"
        source.mkdir()
        output.mkdir()

        files = {
            "Geneformer-V1-10M/model.safetensors": b"synthetic-safetensors",
            "Geneformer-V1-10M/config.json": json.dumps(
                {
                    "architectures": ["BertForMaskedLM"],
                    "hidden_size": 2,
                    "intermediate_size": 4,
                    "max_position_embeddings": 8,
                    "model_type": "bert",
                    "num_attention_heads": 1,
                    "num_hidden_layers": 1,
                    "pad_token_id": 0,
                    "vocab_size": 10,
                },
                sort_keys=True,
            ).encode("utf-8"),
            "geneformer/gene_dictionaries_30m/ensembl_mapping_dict_gc30M.pkl": pickle.dumps(
                {"ENSG000001.1": "ENSG000001"}, protocol=4
            ),
            "geneformer/gene_dictionaries_30m/gene_median_dictionary_gc30M.pkl": pickle.dumps(
                {"ENSG000001": 1.25}, protocol=4
            ),
            "geneformer/gene_dictionaries_30m/gene_name_id_dict_gc30M.pkl": pickle.dumps(
                {"GENE1": "ENSG000001"}, protocol=4
            ),
            "geneformer/gene_dictionaries_30m/token_dictionary_gc30M.pkl": pickle.dumps(
                {"ENSG000001": 7}, protocol=4
            ),
        }
        for relative, payload in files.items():
            path = source / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(payload)

        first_weight = _artifact(
            source / "Geneformer-V1-10M/model.safetensors",
            "Geneformer-V1-10M/model.safetensors",
        )
        placeholder = {
            "path": "not-selected",
            "sha256": "0" * 64,
            "size_bytes": 1,
        }
        manifest = {
            "schema_version": "masld-bench-upstream-checkpoint-preflight-v1",
            "repository": "ctheodoris/Geneformer",
            "repository_revision": "0" * 40,
            "artifacts": [first_weight, placeholder, placeholder],
            "auxiliary_artifacts": {
                "V1": [
                    _artifact(source / relative, relative)
                    for relative in files
                    if relative != "Geneformer-V1-10M/model.safetensors"
                ]
            },
            "architecture_contract": {
                "geneformer_v1_10m": {
                    "architectures": ["BertForMaskedLM"],
                    "hidden_size": 2,
                    "intermediate_size": 4,
                    "max_position_embeddings": 8,
                    "model_type": "bert",
                    "num_attention_heads": 1,
                    "num_hidden_layers": 1,
                    "pad_token_id": 0,
                    "vocab_size": 10,
                }
            },
        }
        manifest_path = root / "checkpoints.json"
        manifest_path.write_text(
            json.dumps(manifest, sort_keys=True, separators=(",", ":")) + "\n",
            encoding="utf-8",
        )
        return source, output, manifest_path

    def test_hash_first_stage_publishes_only_verified_safe_content(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            source, output, manifest = self._fixture(Path(temporary))
            target = stage_geneformer_bundle(
                model_id="geneformer_v1_10m",
                source_root=source,
                output_root=output,
                manifest_path=manifest,
            )
            verify_frozen_tree(target)
            self.assertTrue(
                (target / "upstream/Geneformer-V1-10M/model.safetensors").is_file()
            )
            sanitized = target / (
                "sanitized/geneformer/gene_dictionaries_30m/"
                "token_dictionary_gc30M.pkl"
            )
            self.assertEqual(
                load_restricted_plain_mapping(sanitized.read_bytes()),
                {"ENSG000001": 7},
            )
            self.assertTrue(sanitized.with_suffix(".pkl.json").is_file())
            self.assertFalse(
                (
                    target
                    / "upstream/geneformer/gene_dictionaries_30m/"
                    "token_dictionary_gc30M.pkl"
                ).exists()
            )
            bundle = json.loads((target / "bundle.json").read_text(encoding="utf-8"))
            self.assertIs(bundle["weight_content_executed"], False)
            self.assertIs(bundle["downloads_performed"], False)
            with self.assertRaisesRegex(
                CheckpointPreflightError, "refusing to replace"
            ):
                stage_geneformer_bundle(
                    model_id="geneformer_v1_10m",
                    source_root=source,
                    output_root=output,
                    manifest_path=manifest,
                )

    def test_hash_mismatch_fails_before_publication(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            source, output, manifest = self._fixture(Path(temporary))
            weight = source / "Geneformer-V1-10M/model.safetensors"
            weight.write_bytes(b"tampered")
            with self.assertRaisesRegex(CheckpointPreflightError, "size differs"):
                stage_geneformer_bundle(
                    model_id="geneformer_v1_10m",
                    source_root=source,
                    output_root=output,
                    manifest_path=manifest,
                )
            self.assertEqual(list(output.iterdir()), [])


if __name__ == "__main__":
    unittest.main()
