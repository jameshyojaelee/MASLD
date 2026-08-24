#!/usr/bin/env python3
"""Strictly load and probe the restricted CREsted Enformer port on fixed DNA."""

from __future__ import annotations

import argparse
import csv
import gzip
from hashlib import sha256
import importlib.util
import json
from pathlib import Path, PurePosixPath
import sys
import tarfile
from typing import Any, BinaryIO


OUTER_SHA256 = "628f67f540304d4d0e143176dc824ed72b3413f78f2fa2efe5d4f0ab51ea1bcc"
KERAS_NAME = "enformer_crested_human.keras"
KERAS_SHA256 = "29dc3835c1d13a6c6bd93315b554075107cc409a9179b11ee0e00180f67cef56"
KERAS_SIZE = 985_569_856
ATTENTION_MEMBER = (
    "CREsted-934eea18c97992cd098d547cc09edb20b79c28fd/"
    "src/crested/tl/zoo/utils/_attention.py"
)
ATTENTION_SHA256 = "6c195deb31303d5c1664341308f84625c3d862038fafc8725c781ce33d1c0dac"
LAYERS_MEMBER = (
    "CREsted-934eea18c97992cd098d547cc09edb20b79c28fd/"
    "src/crested/tl/zoo/utils/_layers.py"
)
LAYERS_SHA256 = "44d4befd0e7faca5bf4eea57e94c296147d51818b805d6a8a8a2189ef2f4b457"
LAYERS_SIZE = 23_651
GELU_ENF_DEFINITION = (
    '@keras.saving.register_keras_serializable(package="crested", name="gelu_enf")\n'
    "def gelu_enf(x):\n"
    '    """Very simple gelu approximation, used in Enformer, so needed to get equivalent results."""\n'
    "    return keras.ops.sigmoid(1.702 * x) * x\n"
)
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


class CheckpointProbeError(ValueError):
    """Raised when the restricted checkpoint or fixture differs."""


def _digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def _copy_and_hash(source: BinaryIO, destination: BinaryIO) -> tuple[str, int]:
    value = sha256()
    size = 0
    for block in iter(lambda: source.read(1024 * 1024), b""):
        value.update(block)
        destination.write(block)
        size += len(block)
    return value.hexdigest(), size


def _safe_member(value: str) -> str:
    path = PurePosixPath(value)
    if (
        not value
        or value.startswith("/")
        or "\\" in value
        or "\x00" in value
        or any(part in {"", ".", ".."} for part in path.parts)
    ):
        raise CheckpointProbeError("archive member path is unsafe")
    return path.as_posix()


def _extract_exact_member(
    archive: Path,
    member_name: str,
    destination: Path,
    expected_sha256: str,
    expected_size: int | None = None,
) -> None:
    found = False
    with tarfile.open(archive, mode="r|gz") as handle:
        for member in handle:
            name = _safe_member(member.name)
            if name != member_name:
                continue
            if found or not member.isfile() or member.issym() or member.islnk():
                raise CheckpointProbeError("archive member type or multiplicity differs")
            source = handle.extractfile(member)
            if source is None:
                raise CheckpointProbeError("archive member is unreadable")
            with destination.open("xb") as target:
                observed_sha256, observed_size = _copy_and_hash(source, target)
            if observed_sha256 != expected_sha256 or (
                expected_size is not None and observed_size != expected_size
            ):
                raise CheckpointProbeError("archive member identity differs")
            found = True
    if not found:
        raise CheckpointProbeError("archive member is absent")


def _read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name: str | None = None
    pieces: list[str] = []
    with gzip.open(path, "rt", encoding="ascii") as handle:
        for line in handle:
            value = line.rstrip("\n")
            if value.startswith(">"):
                if name is not None:
                    records.append((name, "".join(pieces).upper()))
                name, pieces = value[1:], []
            elif name is None:
                raise CheckpointProbeError("FASTA sequence precedes its header")
            else:
                pieces.append(value)
    if name is not None:
        records.append((name, "".join(pieces).upper()))
    if not records or len({name for name, _ in records}) != len(records):
        raise CheckpointProbeError("FASTA record identity differs")
    return records


def _one_hot(sequence: str) -> Any:
    import numpy as np

    lookup = np.zeros((256, 4), dtype=np.float32)
    for index, base in enumerate("ACGT"):
        lookup[ord(base), index] = 1.0
    encoded = sequence.encode("ascii")
    array = lookup[np.frombuffer(encoded, dtype=np.uint8)]
    if array.shape != (len(sequence), 4) or not np.isfinite(array).all():
        raise CheckpointProbeError("one-hot conversion differs")
    return array


def _load_fixture_row(manifest: Path, fixture_id: str) -> dict[str, str]:
    with manifest.open("r", encoding="utf-8", newline="") as handle:
        rows = [
            dict(row)
            for row in csv.DictReader(handle, delimiter="\t")
            if row.get("model_id") == "enformer" and row.get("fixture_id") == fixture_id
        ]
    if len(rows) != 1:
        raise CheckpointProbeError("fixture manifest selection differs")
    return rows[0]


def _validate_gelu_source(source: str) -> None:
    if GELU_ENF_DEFINITION not in source:
        raise CheckpointProbeError("CREsted gelu_enf source definition differs")


def _import_custom_layers(attention_path: Path, layers_path: Path) -> dict[str, Any]:
    import keras

    name = "crested.tl.zoo.utils._attention"
    spec = importlib.util.spec_from_file_location(name, attention_path)
    if spec is None or spec.loader is None:
        raise CheckpointProbeError("custom-layer source cannot be imported")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    previous_dont_write_bytecode = sys.dont_write_bytecode
    sys.dont_write_bytecode = True
    try:
        spec.loader.exec_module(module)
    finally:
        sys.dont_write_bytecode = previous_dont_write_bytecode
    source = layers_path.read_text(encoding="utf-8")
    _validate_gelu_source(source)

    @keras.saving.register_keras_serializable(package="crested", name="gelu_enf")
    def gelu_enf(x: Any) -> Any:
        return keras.ops.sigmoid(1.702 * x) * x

    return {
        "AttentionPool1D": module.AttentionPool1D,
        "MultiheadAttention": module.MultiheadAttention,
        "crested>AttentionPool1D": module.AttentionPool1D,
        "crested>MultiheadAttention": module.MultiheadAttention,
        "gelu_enf": gelu_enf,
        "crested>gelu_enf": gelu_enf,
    }


def probe(
    *,
    model_archive: Path,
    source_archive: Path,
    fixture_fasta: Path,
    fixture_manifest: Path,
    fixture_id: str,
    scratch: Path,
    output: Path,
) -> dict[str, Any]:
    import keras
    import numpy as np
    import tensorflow as tf

    if output.exists() or not scratch.is_dir():
        raise CheckpointProbeError("output or scratch contract is invalid")
    if any(
        path.is_symlink() or not path.is_file()
        for path in (model_archive, source_archive, fixture_fasta, fixture_manifest)
    ):
        raise CheckpointProbeError("probe input is not a regular file")
    if _digest(model_archive) != OUTER_SHA256:
        raise CheckpointProbeError("outer model archive checksum differs")

    model_path = scratch / KERAS_NAME
    attention_path = scratch / "crested_attention.py"
    layers_path = scratch / "crested_layers.py"
    _extract_exact_member(
        model_archive, KERAS_NAME, model_path, KERAS_SHA256, KERAS_SIZE
    )
    _extract_exact_member(
        source_archive, ATTENTION_MEMBER, attention_path, ATTENTION_SHA256
    )
    _extract_exact_member(
        source_archive, LAYERS_MEMBER, layers_path, LAYERS_SHA256, LAYERS_SIZE
    )
    custom_objects = _import_custom_layers(attention_path, layers_path)
    tf.keras.utils.set_random_seed(20260824)
    model = keras.models.load_model(
        model_path,
        custom_objects=custom_objects,
        compile=False,
        safe_mode=True,
    )
    model.trainable = False
    if list(model.input_shape) != [None, 196608, 4]:
        raise CheckpointProbeError(f"model input shape differs: {model.input_shape!r}")
    if list(model.output_shape) != [None, 896, 5313]:
        raise CheckpointProbeError(f"model output shape differs: {model.output_shape!r}")

    row = _load_fixture_row(fixture_manifest, fixture_id)
    records = dict(_read_fasta(fixture_fasta))
    sequence = records.get(f"{fixture_id}|enformer|REF")
    if sequence is None or len(sequence) != 196608:
        raise CheckpointProbeError("fixture FASTA sequence differs")
    if sha256(sequence.encode("ascii")).hexdigest() != row["reference_sequence_sha256"]:
        raise CheckpointProbeError("fixture reference sequence checksum differs")
    input_start = int(row["input_start"])
    variant_position0 = int(row["variant_pos_1based"]) - 1 - input_start
    if variant_position0 != 196608 // 2 or sequence[variant_position0] != row["ref"]:
        raise CheckpointProbeError("fixture allele anchoring differs")
    alternative = (
        sequence[:variant_position0]
        + row["alt"]
        + sequence[variant_position0 + 1 :]
    )
    reverse_complement = sequence.translate(COMPLEMENT)[::-1]
    if sha256(alternative.encode("ascii")).hexdigest() != row["alternative_sequence_sha256"]:
        raise CheckpointProbeError("fixture alternative sequence checksum differs")
    if sha256(reverse_complement.encode("ascii")).hexdigest() != row["reverse_complement_reference_sha256"]:
        raise CheckpointProbeError("fixture reverse-complement checksum differs")

    def predict(value: str) -> Any:
        tensor = tf.convert_to_tensor(_one_hot(value)[None, :, :], dtype=tf.float32)
        result = np.asarray(model(tensor, training=False).numpy(), dtype=np.float32)
        if result.shape != (1, 896, 5313) or not np.isfinite(result).all():
            raise CheckpointProbeError("model prediction shape or finiteness differs")
        if np.any(result < 0):
            raise CheckpointProbeError("softplus model emitted a negative value")
        return result[0]

    reference = predict(sequence)
    repeat = predict(sequence)
    if not np.array_equal(reference, repeat):
        raise CheckpointProbeError("deterministic repeat differs")
    alternative_prediction = predict(alternative)
    rc_raw = predict(reverse_complement)
    rc_restored = rc_raw[::-1, :]
    ref_sum = reference.sum(axis=0, dtype=np.float64)
    alt_sum = alternative_prediction.sum(axis=0, dtype=np.float64)
    sad = alt_sum - ref_sum
    sar = np.log2(alt_sum + 1.0) - np.log2(ref_sum + 1.0)
    if not np.isfinite(sad).all() or not np.isfinite(sar).all():
        raise CheckpointProbeError("variant summary is non-finite")

    output.mkdir(parents=True, mode=0o750)
    predictions_path = output / "native_fixture_predictions.npz"
    np.savez_compressed(
        predictions_path,
        reference=reference,
        alternative=alternative_prediction,
        reverse_complement_raw=rc_raw,
        reverse_complement_restored=rc_restored,
        sad_alt_minus_ref=sad,
        sar_alt_minus_ref=sar,
    )
    with np.load(predictions_path, allow_pickle=False) as frozen:
        if set(frozen.files) != {
            "reference",
            "alternative",
            "reverse_complement_raw",
            "reverse_complement_restored",
            "sad_alt_minus_ref",
            "sar_alt_minus_ref",
        }:
            raise CheckpointProbeError("frozen prediction keys differ")

    devices = tf.config.list_physical_devices("GPU")
    receipt = {
        "schema_version": "masld-bench-enformer-crested-checkpoint-probe-v1",
        "status": "pass",
        "model_id": "enformer_crested_restricted_port",
        "registered_scientific_identity": "restricted_conversion_not_native_sonnet",
        "keras_safe_mode": True,
        "compile": False,
        "custom_layer_source_sha256": ATTENTION_SHA256,
        "custom_activation_source_sha256": LAYERS_SHA256,
        "checkpoint_member_sha256": KERAS_SHA256,
        "tensorflow": tf.__version__,
        "keras": keras.__version__,
        "gpu_devices": [device.name for device in devices],
        "parameter_count": int(model.count_params()),
        "input_shape": list(model.input_shape),
        "output_shape": list(model.output_shape),
        "fixture_id": fixture_id,
        "fixture_genomic_fold": int(row["genomic_fold"]),
        "fixture_reference_sequence_sha256": row["reference_sequence_sha256"],
        "fixture_alternative_sequence_sha256": row["alternative_sequence_sha256"],
        "allele_delta_sign": "ALT_minus_REF",
        "deterministic_repeat_bit_identical": True,
        "reference_prediction_sha256": sha256(reference.tobytes()).hexdigest(),
        "alternative_prediction_sha256": sha256(alternative_prediction.tobytes()).hexdigest(),
        "reverse_complement_raw_prediction_sha256": sha256(rc_raw.tobytes()).hexdigest(),
        "sad_sha256": sha256(sad.tobytes()).hexdigest(),
        "sar_sha256": sha256(sar.tobytes()).hexdigest(),
        "predictions_artifact_sha256": _digest(predictions_path),
        "all_outputs_finite_and_nonnegative": True,
        "project_data_read": False,
        "outcomes_read": False,
        "sealed_features_or_labels_read": False,
        "native_sonnet_parity_established": False,
        "training_or_adaptation_executed": False,
        "full_fine_tuning_allowed_under_current_terms": False,
        "open_champion_eligible": False,
        "allowed_use": "internal_restricted_static_track_and_variant_component_comparator",
        "terminal_disposition": "restricted_port_strict_restore_and_fixed_forward_passed_native_sonnet_parity_unresolved",
    }
    (output / "receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    model_path.unlink()
    attention_path.unlink()
    layers_path.unlink()
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model-archive", type=Path, required=True)
    parser.add_argument("--source-archive", type=Path, required=True)
    parser.add_argument("--fixture-fasta", type=Path, required=True)
    parser.add_argument("--fixture-manifest", type=Path, required=True)
    parser.add_argument("--fixture-id", required=True)
    parser.add_argument("--scratch", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    result = probe(
        model_archive=arguments.model_archive,
        source_archive=arguments.source_archive,
        fixture_fasta=arguments.fixture_fasta,
        fixture_manifest=arguments.fixture_manifest,
        fixture_id=arguments.fixture_id,
        scratch=arguments.scratch,
        output=arguments.output,
    )
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
