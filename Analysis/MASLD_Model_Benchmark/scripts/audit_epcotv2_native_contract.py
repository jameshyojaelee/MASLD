#!/usr/bin/env python3
"""Audit EPCOTv2's source-defined observed-ATAC requirements without loading weights."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import tarfile
from typing import Any


class EPCOTv2AuditError(RuntimeError):
    """Raised when frozen EPCOTv2 evidence differs from the included requirements."""


CODE_REVISION = "f20362ab5139ac53448745cdca3480f952db5cf1"
SPACE_REVISION = "f93a2a944a49797a111f0294709416d38bedbcbe"
REFERENCE_REVISION = "f445058917800811cb82c38b2c806674cc02b47f"
CHECKPOINT_SHA256 = "b824a80ed238e64e15c5a6eeae91329d24f04825c62d2a9dd328658822aba45d"
ADMISSION_ARTIFACTS_SHA256 = "26e12915b05c2727c45eca091f85ec89f96c205e817c6ec9b64f181a563e7577"
LEGACY_FIXTURE_ARTIFACTS_SHA256 = "4ac1e23e884dab0318c2c2fbf18e5a8ea244c5db8d0bf40ef95f3479b4f25607"

PUBLICATION_FILES = {
    "LICENSE": "2aff8b165075fc1c015e2b0f9df0adb319e965ab51cf3806c9f99dadec797b7e",
    "README.md": "738e2190e60e5217a6076c5b014114025c280685946c746127de7a19ba89fd05",
    "data/README.md": "aebaa4fedb52d138a3864963bd16a595f47daf08e1327176399aa06bb31b64b2",
    "data/input_region_600kb.bed": "d6991d59eb7469d014be9583e15d6d2799ef3af83154e4eea7648932193605d7",
    "src/atac_process.py": "7f1680e0076a59ee98a7b11d308d2f22a687d6614d9c4677a24b3bafa13ec272",
    "src/README.md": "3fca92e9e1b088dcb46621e488fe9a592b3fbd8a6148546dcfe27c9e058ab875",
    "src/model.py": "853ebbb51c4f68dbf0d6189218a18aec6961afdc83957973b049958efb37143a",
    "src/train.py": "b44ff8c6a0db1873ca339454e51fdbac80f70b199f6f4a3853ef59ff7aaeb40d",
    "src/util.py": "88bfc283336df4d79d6fd897fd20466203d192ba4db3a096a94bd698e88b251b",
}

SPACE_FILES = {
    "epcotv2_space_curriculum_lora_prompt_model.py": {
        "transport_sha256": "2445cc49ce1e71eff25da2722df3e2812e0c75ead069d10ae784d95670804640",
        "git_blob_content_sha256": "afdaeef02743cfb9626f7a8deaef39daf4a13679f3e182744d8fd51444c6b7d6",
    },
    "epcotv2_space_func_gradio.py": {
        "transport_sha256": "b5f53e0811e2e270e40b6506d863a32b16252828e506c8cf98dbec20d6a825dc",
        "git_blob_content_sha256": "529796f944c49d8aa056cd06d249c24d7449da938386daebea545936f4adc419",
    },
    "epcotv2_space_erna_input_region_600kb.bed": {
        "transport_sha256": "d6991d59eb7469d014be9583e15d6d2799ef3af83154e4eea7648932193605d7",
        "git_blob_content_sha256": "302557ff2cbc6cd09829afc1dc97e702e42818645ac533070f51c1bba1fcd103",
    },
}


def digest_bytes(value: bytes) -> str:
    return sha256(value).hexdigest()


def digest_path(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise EPCOTv2AuditError(f"expected JSON object: {path}")
    return value


def verify_artifact(root: Path, expected_manifest_sha256: str) -> None:
    from masld_bench.artifacts import ArtifactError, verify_frozen_tree

    try:
        verify_frozen_tree(root)
    except ArtifactError as error:
        raise EPCOTv2AuditError(f"frozen artifact differs: {root}: {error}") from error
    if digest_path(root / "ARTIFACTS.json") != expected_manifest_sha256:
        raise EPCOTv2AuditError(f"unexpected frozen artifact identity: {root}")


def read_tar_member(archive: Path, relative: str) -> bytes:
    member_name = f"general_AI_model-{CODE_REVISION}/{relative}"
    with tarfile.open(archive, "r:gz") as handle:
        try:
            member = handle.getmember(member_name)
        except KeyError as error:
            raise EPCOTv2AuditError(f"missing publication source: {relative}") from error
        if not member.isfile():
            raise EPCOTv2AuditError(f"publication source is not a regular file: {relative}")
        stream = handle.extractfile(member)
        if stream is None:
            raise EPCOTv2AuditError(f"cannot read publication source: {relative}")
        return stream.read()


def api_terms(metadata: dict[str, Any], expected_revision: str | None = None) -> dict[str, Any]:
    revision = metadata.get("sha")
    if expected_revision is not None and revision != expected_revision:
        raise EPCOTv2AuditError(
            f"Hugging Face revision differs: {revision!r} != {expected_revision!r}"
        )
    card = metadata.get("cardData") or {}
    siblings = metadata.get("siblings") or []
    if not isinstance(card, dict) or not isinstance(siblings, list):
        raise EPCOTv2AuditError("Hugging Face metadata structure differs")
    license_files = sorted(
        str(row.get("rfilename"))
        for row in siblings
        if isinstance(row, dict)
        and str(row.get("rfilename", "")).lower()
        in {"license", "license.md", "license.txt", "copying"}
    )
    declared_license = card.get("license")
    return {
        "revision": revision,
        "declared_license": declared_license,
        "license_files": license_files,
        "terms_declared": bool(declared_license or license_files),
    }


def pad_sequence(matrix: Any, pad_len: int = 300) -> Any:
    """Reproduce the MIT source's neighbor-bin padding on anonymous basis channels."""
    import numpy as np

    if matrix.ndim != 3 or matrix.shape[1:] != (4, 1000):
        raise EPCOTv2AuditError("sequence basis matrix must be bins_by_4_by_1000")
    paddings = np.zeros((1, 4, pad_len), dtype=np.int8)
    downstream = np.concatenate((paddings, matrix[:, :, -pad_len:]), axis=0)[:-1]
    upstream = np.concatenate((matrix[:, :, :pad_len], paddings), axis=0)[1:]
    return np.concatenate((downstream, matrix, upstream), axis=2)


def pad_signal(matrix: Any, pad_len: int = 300) -> Any:
    """Reproduce the MIT source's neighbor-bin padding for observed ATAC."""
    import numpy as np

    if matrix.ndim != 2 or matrix.shape[1] != 1000:
        raise EPCOTv2AuditError("ATAC matrix must be bins_by_1000")
    paddings = np.zeros(pad_len, dtype=np.float32)
    downstream = np.vstack((paddings, matrix[:, -pad_len:]))[:-1]
    upstream = np.vstack((matrix[:, :pad_len], paddings))[1:]
    return np.hstack((downstream, matrix, upstream))


def build_source_defined_fixture(output: Path) -> dict[str, Any]:
    import numpy as np

    bins = 602
    positions = np.arange(bins * 1000, dtype=np.int64).reshape(bins, 1000)
    basis = np.zeros((bins, 4, 1000), dtype=np.int8)
    for channel in range(4):
        basis[:, channel, :] = positions % 4 == channel
    raw_atac = ((positions % 251) + 1).astype(np.float32) / 252.0

    raw_locus_mask = np.zeros_like(raw_atac, dtype=np.bool_)
    raw_locus_mask[301, :] = True
    masked_raw_atac = raw_atac.copy()
    masked_raw_atac[raw_locus_mask] = 0.0

    padded_basis = pad_sequence(basis)[1:601]
    padded_atac = pad_signal(masked_raw_atac)[1:601, None, :]
    propagated_mask = pad_signal(raw_locus_mask.astype(np.float32))[1:601, None, :] > 0
    tensor = np.concatenate((padded_basis.astype(np.float32), padded_atac), axis=1)

    if tensor.shape != (600, 5, 1600) or propagated_mask.shape != (600, 1, 1600):
        raise EPCOTv2AuditError("source-defined fixture shape differs")
    if not np.all(tensor[:, :4].sum(axis=1) == 1):
        raise EPCOTv2AuditError("anonymous sequence basis is not one-hot")
    if int(propagated_mask.sum()) != 1600 or not np.all(tensor[:, 4:5][propagated_mask] == 0):
        raise EPCOTv2AuditError("raw-locus mask did not propagate through padded ATAC")
    affected_bins = np.flatnonzero(propagated_mask.any(axis=(1, 2))).tolist()
    if affected_bins != [299, 300, 301]:
        raise EPCOTv2AuditError("raw-locus mask affected unexpected model bins")

    path = output / "native_preprocessing_fixture.npz"
    np.savez_compressed(
        path,
        sequence_plus_observed_atac=tensor[None],
        propagated_raw_atac_locus_mask=propagated_mask[None],
    )
    return {
        "path": path.name,
        "sha256": digest_path(path),
        "tensor_shape": [1, 600, 5, 1600],
        "tensor_dtype": "float32",
        "sequence_channel_semantics": ["basis_0", "basis_1", "basis_2", "basis_3"],
        "base_to_channel_mapping_claimed": False,
        "atac_raw_locus_mask_bp": 1000,
        "atac_mask_propagated_model_bins": affected_bins,
        "atac_mask_occurrences_after_padding": int(propagated_mask.sum()),
        "model_checkpoint_loaded": False,
        "model_forward_executed": False,
    }


def require_source_fragments(sources: dict[str, bytes]) -> None:
    requirements = {
        "src/README.md": [
            b"--normalizeUsing RPGC",
            b"--effectiveGenomeSize 2913022398",
            b"--Offset 1",
            b"--binSize 1",
        ],
        "src/util.py": [
            b"def pad_signal_matrix(matrix, pad_len=300):",
            b"def pad_seq_matrix(matrix, pad_len=300):",
            b"reshape(4, -1, 1000).swapaxes(0, 1)",
            b"reshape(-1, 1000)",
        ],
        "src/train.py": [
            b"def split_dataset(seed=24):",
            b"valid_split = int(np.floor(dataset_size * 0.8))",
            b"test_split = int(np.floor(dataset_size * 0.9))",
            b"train_indices, valid_indices= indices[:valid_split], indices[test_split:]",
        ],
        "src/model.py": [
            b"x_seq,x_atac=x[:,:4,:],x[:,4:,:]",
            b"x_rep_update[:, self.crop:-self.crop, :]",
            b"nn.Linear(256,247)",
            b"nn.Linear(720,708)",
        ],
    }
    for name, fragments in requirements.items():
        for fragment in fragments:
            if fragment not in sources[name]:
                raise EPCOTv2AuditError(
                    f"source-defined contract fragment differs: {name}: {fragment!r}"
                )


def audit(args: argparse.Namespace) -> dict[str, Any]:
    verify_artifact(args.admission, ADMISSION_ARTIFACTS_SHA256)
    verify_artifact(args.legacy_fixture, LEGACY_FIXTURE_ARTIFACTS_SHA256)
    receipt = load_json(args.admission / "admission_sources_receipt.json")
    if receipt["epcotv2"]["checkpoint_lfs_sha256"] != CHECKPOINT_SHA256:
        raise EPCOTv2AuditError("checkpoint identity differs")
    if receipt["epcotv2"]["checkpoint_object_downloaded"]:
        raise EPCOTv2AuditError("checkpoint object crossed the terms gate")

    archive = args.admission / "downloads/epcotv2_publication_code.tar.gz"
    sources = {name: read_tar_member(archive, name) for name in PUBLICATION_FILES}
    for name, expected in PUBLICATION_FILES.items():
        if digest_bytes(sources[name]) != expected:
            raise EPCOTv2AuditError(f"publication source identity differs: {name}")
    require_source_fragments(sources)

    metadata_root = args.admission / "metadata"
    space_records = {}
    for name, expected in SPACE_FILES.items():
        content = (metadata_root / name).read_bytes()
        transport_sha = digest_bytes(content)
        logical_sha = digest_bytes(content[:-1] if content.endswith(b"\n") else content)
        if transport_sha != expected["transport_sha256"] or logical_sha != expected["git_blob_content_sha256"]:
            raise EPCOTv2AuditError(f"Space source identity differs: {name}")
        space_records[name] = {
            "transport_sha256": transport_sha,
            "git_blob_content_sha256": logical_sha,
            "transport_appended_terminal_lf": content.endswith(b"\n"),
        }

    publication_model = sources["src/model.py"]
    space_model = (metadata_root / "epcotv2_space_curriculum_lora_prompt_model.py").read_bytes()
    publication_active_lines = [
        line.strip()
        for line in publication_model.splitlines()
        if not line.lstrip().startswith(b"#")
    ]
    space_active_lines = [
        line.strip()
        for line in space_model.splitlines()
        if not line.lstrip().startswith(b"#")
    ]
    dropout_fragment = b"F.dropout(x_seq_1,0.05,training=self.training)"
    source_divergence = {
        "byte_identical": publication_model == space_model,
        "publication_training_dropout_0_05_present": any(
            dropout_fragment in line for line in publication_active_lines
        ),
        "space_training_dropout_0_05_active": any(
            dropout_fragment in line for line in space_active_lines
        ),
        "space_prompt_path_present": b"self.prompt_token_weight" in space_model,
        "space_optional_external_heads_present": b"self.external=external" in space_model,
    }
    if source_divergence != {
        "byte_identical": False,
        "publication_training_dropout_0_05_present": True,
        "space_training_dropout_0_05_active": False,
        "space_prompt_path_present": True,
        "space_optional_external_heads_present": True,
    }:
        raise EPCOTv2AuditError("publication-to-Space source divergence differs")

    current_space = api_terms(load_json(args.terms_snapshot / "space_current.json"))
    pinned_space = api_terms(
        load_json(args.terms_snapshot / "space_pinned.json"), SPACE_REVISION
    )
    current_reference = api_terms(load_json(args.terms_snapshot / "reference_current.json"))
    pinned_reference = api_terms(
        load_json(args.terms_snapshot / "reference_pinned.json"), REFERENCE_REVISION
    )
    endpoint_status = load_json(args.terms_snapshot / "license_endpoint_status.json")
    if pinned_space["terms_declared"] or pinned_reference["terms_declared"]:
        raise EPCOTv2AuditError("pinned terms status changed; requires manual reassessment")

    legacy = load_json(args.legacy_fixture / "fixtures/epcotv2/receipt.json")
    if legacy["model_checkpoint_loaded"] or legacy["model_forward_executed"]:
        raise EPCOTv2AuditError("legacy fixture crossed the checkpoint terms gate")

    args.output.mkdir(parents=True, exist_ok=False)
    fixture = build_source_defined_fixture(args.output)
    result = {
        "schema_version": "masld-bench-epcotv2-native-contract-audit-v1",
        "status": "pass",
        "checkpoint": {
            "sha256": CHECKPOINT_SHA256,
            "object_downloaded": False,
            "deserialized": False,
        },
        "source": {
            "publication_code_revision": CODE_REVISION,
            "publication_code_license": "MIT",
            "publication_files_sha256": PUBLICATION_FILES,
            "space_revision": SPACE_REVISION,
            "space_files": space_records,
            "publication_to_space_divergence": source_divergence,
        },
        "terms": {
            "space_current": current_space,
            "space_pinned": pinned_space,
            "reference_current": current_reference,
            "reference_pinned": pinned_reference,
            "license_endpoint_http_status": endpoint_status,
            "checkpoint_weight_terms": "UNDECLARED",
            "terminal_disposition": "blocked_terms_no_checkpoint_download_or_execution",
        },
        "native_input": {
            "model_tensor_shape": ["batch", 600, 5, 1600],
            "sequence_channels": 4,
            "observed_atac_channels": 1,
            "sequence_base_to_channel_mapping": "UNRESOLVED",
            "sequence_reference_identity": "UNRESOLVED",
            "atac_processing": {
                "bam_coverage_normalization": "RPGC",
                "effective_genome_size": 2913022398,
                "offset": 1,
                "source_bin_bp": 1,
                "model_bin_bp": 1000,
                "neighbor_padding_bp_each_side": 300,
            },
        },
        "native_output": {
            "central_output_bins": 500,
            "central_output_span_bp": 500000,
            "modality_groups": 14,
            "atac_is_output": False,
            "same_nucleus_rna_atac_task": "observed_ATAC_plus_sequence_to_RNA_profile",
            "biological_replication_unit": "donor",
            "cell_level_testing_permitted": False,
        },
        "split_reconstruction": {
            "seed": 24,
            "training_fraction": 0.8,
            "validation_fraction": 0.1,
            "unused_middle_fraction": 0.1,
            "checkpoint_binding": "UNRESOLVED",
        },
        "fixture": fixture,
        "fixture_scope": {
            "source_defined_preprocessing_fixture_without_biological_channel_names": "admissible",
            "biologically_interpretable_sequence_fixture": "blocked_reference_identity_and_channel_mapping",
            "checkpoint_forward_or_numeric_fixture": "blocked_exact_code_weight_and_reference_terms",
            "legacy_tensor_local_variant_mask": "insufficient_for_raw_locus_mask_due_to_neighbor_padding",
            "required_variant_mask_rule": "mask_the_raw_genomic_ATAC_locus_before_300bp_neighbor_padding",
        },
        "eligibility": {
            "observed_multiome_development_lane": "preprocessing_fixture_ready_checkpoint_execution_blocked_terms",
            "rna_conditioned_atac_lane": "not_applicable",
            "sealed_champion": "ineligible_without_topology_matched_external_ATAC_and_terms_clearance",
            "gpu_queue_item": False,
        },
        "project_data_read": False,
        "outcomes_read": False,
        "sealed_labels_read": False,
    }
    (args.output / "audit.json").write_text(
        json.dumps(result, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--admission", type=Path, required=True)
    parser.add_argument("--legacy-fixture", type=Path, required=True)
    parser.add_argument("--terms-snapshot", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise EPCOTv2AuditError("output already exists")
    result = audit(args)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
