#!/usr/bin/env python3
"""Resolve EPCOTv2 inclusion evidence without loading weights or biological data."""

from __future__ import annotations

import argparse
import ast
from hashlib import sha256
import json
from pathlib import Path
import tarfile
from typing import Any


class EPCOTv2ResolutionError(RuntimeError):
    """Raised when a released source differs from the fail-closed requirement."""


CODE_REVISION = "f20362ab5139ac53448745cdca3480f952db5cf1"
SPACE_REVISION = "f93a2a944a49797a111f0294709416d38bedbcbe"
REFERENCE_REVISION = "f445058917800811cb82c38b2c806674cc02b47f"
CHECKPOINT_LFS_SHA256 = (
    "b824a80ed238e64e15c5a6eeae91329d24f04825c62d2a9dd328658822aba45d"
)
CHECKPOINT_SIZE_BYTES = 468_200_299
ADMISSION_ARTIFACTS_SHA256 = (
    "26e12915b05c2727c45eca091f85ec89f96c205e817c6ec9b64f181a563e7577"
)
KIPOISEQ_VERSION = "0.5.2"
KIPOISEQ_SDIST_SHA256 = (
    "b1d8f717326522e4fd88baca8bb288050a6da2ec75e23b82006c55c6bc68204d"
)
KIPOISEQ_SDIST_SIZE_BYTES = 34_170
UCSC_HG38_FASTA_URL = (
    "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
)

PUBLICATION_FILES = {
    "LICENSE": "2aff8b165075fc1c015e2b0f9df0adb319e965ab51cf3806c9f99dadec797b7e",
    "requirements.txt": "2bb2625debe9305e5a88ff7decd54e6c57c237e858191834b4db81d06c16dc1c",
    "src/tutorial_utils.py": "afd551fa4d183e2f61a8bfd59403460ff4d55ef26dbabbf723ad178c9287257f",
    "src/util.py": "88bfc283336df4d79d6fd897fd20466203d192ba4db3a096a94bd698e88b251b",
    "epcotv2_basic_tutorial.ipynb": "372671c4e7561a36c19650579af857eb2f879b9ea1f870b0c657581f86c5da2e",
}

SPACE_SOURCE_FILES = (
    "README.md",
    "app.py",
    "curriculum/lora_prompt_model.py",
    "curriculum/loralib/__init__.py",
    "curriculum/loralib/layers.py",
    "curriculum/loralib/utils.py",
    "erna/util.py",
    "func_gradio.py",
    "requirements.txt",
)


def digest_bytes(value: bytes) -> str:
    return sha256(value).hexdigest()


def digest_path(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise EPCOTv2ResolutionError(f"expected JSON object: {path}")
    return value


def verify_frozen_admission(root: Path) -> None:
    from masld_bench.artifacts import ArtifactError, verify_frozen_tree

    try:
        verify_frozen_tree(root)
    except ArtifactError as error:
        raise EPCOTv2ResolutionError(
            f"frozen admission artifact differs: {root}: {error}"
        ) from error
    if digest_path(root / "ARTIFACTS.json") != ADMISSION_ARTIFACTS_SHA256:
        raise EPCOTv2ResolutionError("unexpected admission artifact identity")


def read_exact_tar_member(archive: Path, relative: str) -> bytes:
    member_name = f"general_AI_model-{CODE_REVISION}/{relative}"
    with tarfile.open(archive, "r:gz") as handle:
        try:
            member = handle.getmember(member_name)
        except KeyError as error:
            raise EPCOTv2ResolutionError(
                f"missing publication source: {relative}"
            ) from error
        if not member.isfile():
            raise EPCOTv2ResolutionError(
                f"publication source is not a regular file: {relative}"
            )
        stream = handle.extractfile(member)
        if stream is None:
            raise EPCOTv2ResolutionError(f"cannot read publication source: {relative}")
        return stream.read()


def read_unique_tar_suffix(archive: Path, suffix: str) -> bytes:
    with tarfile.open(archive, "r:gz") as handle:
        matches = [
            member
            for member in handle.getmembers()
            if member.isfile() and member.name.endswith(suffix)
        ]
        if len(matches) != 1:
            raise EPCOTv2ResolutionError(
                f"expected one sdist member ending {suffix!r}; found {len(matches)}"
            )
        stream = handle.extractfile(matches[0])
        if stream is None:
            raise EPCOTv2ResolutionError(f"cannot read sdist member: {matches[0].name}")
        return stream.read()


def literal_assignment(source: bytes, name: str) -> Any:
    tree = ast.parse(source.decode("utf-8"))
    for node in tree.body:
        if not isinstance(node, (ast.Assign, ast.AnnAssign)):
            continue
        targets = node.targets if isinstance(node, ast.Assign) else [node.target]
        if any(isinstance(target, ast.Name) and target.id == name for target in targets):
            try:
                return ast.literal_eval(node.value)
            except (TypeError, ValueError) as error:
                raise EPCOTv2ResolutionError(
                    f"assignment {name} is not a literal"
                ) from error
    raise EPCOTv2ResolutionError(f"assignment not found: {name}")


def verify_one_hot_dna_contract(source: bytes) -> dict[str, Any]:
    tree = ast.parse(source.decode("utf-8"))
    function = next(
        (
            node
            for node in tree.body
            if isinstance(node, ast.FunctionDef) and node.name == "one_hot_dna"
        ),
        None,
    )
    if function is None:
        raise EPCOTv2ResolutionError("kipoiseq one_hot_dna function is absent")
    calls = [node for node in ast.walk(function) if isinstance(node, ast.Call)]
    call = next(
        (
            node
            for node in calls
            if isinstance(node.func, ast.Name) and node.func.id == "one_hot"
        ),
        None,
    )
    if call is None:
        raise EPCOTv2ResolutionError("one_hot_dna does not call one_hot")
    keywords = {keyword.arg: keyword.value for keyword in call.keywords if keyword.arg}
    alphabet = keywords.get("alphabet")
    if not isinstance(alphabet, ast.Name) or alphabet.id != "DNA":
        raise EPCOTv2ResolutionError("one_hot_dna alphabet is not the DNA constant")
    try:
        neutral_alphabet = ast.literal_eval(keywords["neutral_alphabet"])
        neutral_value = ast.literal_eval(keywords["neutral_value"])
    except (KeyError, TypeError, ValueError) as error:
        raise EPCOTv2ResolutionError(
            "one_hot_dna neutral-base contract is not explicit"
        ) from error
    if neutral_alphabet not in (["N"], ("N",)) or float(neutral_value) != 0.25:
        raise EPCOTv2ResolutionError("one_hot_dna neutral-base contract differs")
    return {"neutral_alphabet": ["N"], "neutral_value_per_channel": 0.25}


def api_terms(metadata: dict[str, Any], expected_revision: str | None) -> dict[str, Any]:
    revision = metadata.get("sha")
    if expected_revision is not None and revision != expected_revision:
        raise EPCOTv2ResolutionError(
            f"Hugging Face revision differs: {revision!r} != {expected_revision!r}"
        )
    card = metadata.get("cardData") or {}
    siblings = metadata.get("siblings") or []
    if not isinstance(card, dict) or not isinstance(siblings, list):
        raise EPCOTv2ResolutionError("Hugging Face metadata structure differs")
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


def require_fragments(source: bytes, fragments: tuple[bytes, ...], name: str) -> None:
    for fragment in fragments:
        if fragment not in source:
            raise EPCOTv2ResolutionError(f"source fragment differs: {name}: {fragment!r}")


def audit(args: argparse.Namespace) -> dict[str, Any]:
    verify_frozen_admission(args.admission)
    admission_receipt = load_json(args.admission / "admission_sources_receipt.json")
    checkpoint = admission_receipt["epcotv2"]
    if checkpoint["checkpoint_lfs_sha256"] != CHECKPOINT_LFS_SHA256:
        raise EPCOTv2ResolutionError("checkpoint LFS identity differs")
    if checkpoint["checkpoint_size_bytes"] != CHECKPOINT_SIZE_BYTES:
        raise EPCOTv2ResolutionError("checkpoint LFS size differs")
    if checkpoint["checkpoint_object_downloaded"]:
        raise EPCOTv2ResolutionError("checkpoint object crossed the terms gate")

    publication_archive = (
        args.admission / "downloads/epcotv2_publication_code.tar.gz"
    )
    publication_sources = {
        name: read_exact_tar_member(publication_archive, name)
        for name in PUBLICATION_FILES
    }
    for name, expected in PUBLICATION_FILES.items():
        if digest_bytes(publication_sources[name]) != expected:
            raise EPCOTv2ResolutionError(f"publication source identity differs: {name}")

    requirements = publication_sources["requirements.txt"].decode("utf-8").splitlines()
    if f"kipoiseq=={KIPOISEQ_VERSION}" not in requirements:
        raise EPCOTv2ResolutionError("publication kipoiseq pin differs")
    tutorial_utils = publication_sources["src/tutorial_utils.py"]
    require_fragments(
        tutorial_utils,
        (
            b"kipoiseq.transforms.functional.one_hot_dna(sequence)",
            b"reshape(-1, 1000, 4).swapaxes(1, 2)",
        ),
        "publication src/tutorial_utils.py",
    )
    tutorial = json.loads(
        publication_sources["epcotv2_basic_tutorial.ipynb"].decode("utf-8")
    )
    tutorial_text = json.dumps(tutorial, sort_keys=True)
    if UCSC_HG38_FASTA_URL not in tutorial_text:
        raise EPCOTv2ResolutionError("publication tutorial reference URL differs")

    pypi = load_json(args.source_cache / "metadata/kipoiseq-0.5.2-pypi.json")
    if pypi.get("info", {}).get("version") != KIPOISEQ_VERSION:
        raise EPCOTv2ResolutionError("PyPI kipoiseq version differs")
    sdist_rows = [
        row
        for row in pypi.get("urls", [])
        if isinstance(row, dict) and row.get("packagetype") == "sdist"
    ]
    if len(sdist_rows) != 1:
        raise EPCOTv2ResolutionError("PyPI kipoiseq sdist cardinality differs")
    sdist_row = sdist_rows[0]
    if (
        sdist_row.get("digests", {}).get("sha256") != KIPOISEQ_SDIST_SHA256
        or sdist_row.get("size") != KIPOISEQ_SDIST_SIZE_BYTES
    ):
        raise EPCOTv2ResolutionError("PyPI kipoiseq sdist metadata differs")
    sdist = args.source_cache / "downloads/kipoiseq-0.5.2.tar.gz"
    if digest_path(sdist) != KIPOISEQ_SDIST_SHA256:
        raise EPCOTv2ResolutionError("downloaded kipoiseq sdist identity differs")
    kipoiseq_utils = read_unique_tar_suffix(sdist, "/kipoiseq/utils.py")
    kipoiseq_functional = read_unique_tar_suffix(
        sdist, "/kipoiseq/transforms/functional.py"
    )
    dna_alphabet = literal_assignment(kipoiseq_utils, "DNA")
    if dna_alphabet != ["A", "C", "G", "T"]:
        raise EPCOTv2ResolutionError(f"kipoiseq DNA alphabet differs: {dna_alphabet!r}")
    neutral_contract = verify_one_hot_dna_contract(kipoiseq_functional)

    space_api = load_json(args.source_cache / "metadata/space-pinned.json")
    reference_api = load_json(args.source_cache / "metadata/reference-pinned.json")
    space_current_api = load_json(args.source_cache / "metadata/space-current.json")
    reference_current_api = load_json(
        args.source_cache / "metadata/reference-current.json"
    )
    space_terms = api_terms(space_api, SPACE_REVISION)
    reference_terms = api_terms(reference_api, REFERENCE_REVISION)
    if space_terms["terms_declared"] or reference_terms["terms_declared"]:
        raise EPCOTv2ResolutionError("pinned EPCOTv2 terms changed; manual review required")

    siblings = {
        str(row.get("rfilename"))
        for row in space_api.get("siblings", [])
        if isinstance(row, dict)
    }
    if not set(SPACE_SOURCE_FILES).issubset(siblings):
        raise EPCOTv2ResolutionError("pinned Space source closure differs")
    space_root = args.source_cache / "sources/space"
    space_sources = {
        name: (space_root / name).read_bytes() for name in SPACE_SOURCE_FILES
    }
    source_hashes = {name: digest_bytes(value) for name, value in space_sources.items()}

    require_fragments(
        space_sources["app.py"],
        (b"from func_gradio import run_epcotv2", b"gr.Interface(fn=run_epcotv2"),
        "Space app.py",
    )
    require_fragments(
        space_sources["func_gradio.py"],
        (
            b'model.load_state_dict(torch.load(os.path.join(sys.path[0], "models/human_model.pt"), map_location=torch.device(\'cpu\')))',
            b"pickle.load(f)",
            b"pickle.dump(pred_outputs, tmp)",
            b"load_ref_genome(chrom)",
        ),
        "Space func_gradio.py",
    )
    require_fragments(
        space_sources["erna/util.py"],
        (
            b'repo_id="luosanj/epcotv2_data"',
            b"hf_hub_download(",
            b"reshape(4, -1, 1000).swapaxes(0, 1)",
        ),
        "Space erna/util.py",
    )
    if b"revision=" in space_sources["erna/util.py"]:
        raise EPCOTv2ResolutionError(
            "Space reference fetch gained a revision argument; manual review required"
        )
    space_requirements = space_sources["requirements.txt"].decode("utf-8").splitlines()
    if "torch==2.2.1" not in space_requirements or "huggingface_hub" not in space_requirements:
        raise EPCOTv2ResolutionError("Space runtime dependency contract differs")

    dataset_readme = (args.source_cache / "sources/reference/README.md").read_text(
        encoding="utf-8"
    )
    if "[More Information Needed]" not in dataset_readme:
        raise EPCOTv2ResolutionError("reference dataset documentation changed")
    reference_siblings = {
        str(row.get("rfilename"))
        for row in reference_api.get("siblings", [])
        if isinstance(row, dict)
    }
    expected_reference_files = {f"chr{chrom}.npz" for chrom in range(1, 23)} | {
        "chrX.npz"
    }
    if not expected_reference_files.issubset(reference_siblings):
        raise EPCOTv2ResolutionError("reference array chromosome set differs")

    license_status = load_json(
        args.source_cache / "metadata/license-endpoint-status.json"
    )
    if any(status != 404 for status in license_status.values()):
        raise EPCOTv2ResolutionError(
            f"license endpoint status changed: {license_status}; manual review required"
        )
    result = {
        "schema_version": "masld-bench-epcotv2-biological-admission-resolution-v1",
        "audit_status": "passed_fail_closed",
        "scientific_lane": {
            "native_input": "DNA_sequence_plus_observed_ATAC",
            "native_output": "multi_assay_profiles_and_contacts_including_RNA",
            "rna_conditioned_atac": "not_applicable",
            "observed_multiome_lane": "applicable_after_all_admission_gates_pass",
        },
        "publication_tutorial_sequence_contract": {
            "code_revision": CODE_REVISION,
            "code_license": "MIT",
            "kipoiseq_version": KIPOISEQ_VERSION,
            "kipoiseq_sdist_sha256": KIPOISEQ_SDIST_SHA256,
            "base_to_channel": {base: index for index, base in enumerate(dna_alphabet)},
            "neutral_base": neutral_contract,
            "channel_mapping_status": "resolved_for_publication_tutorial_path",
        },
        "space_reference_array_contract": {
            "space_revision": SPACE_REVISION,
            "reference_repository_revision_inspected": REFERENCE_REVISION,
            "reference_files": sorted(expected_reference_files),
            "stored_shape_transform": "sparse_rows_to_4_by_genome_then_bins_by_4_by_1000",
            "base_to_row_mapping": "UNRESOLVED",
            "fasta_source_and_checksum": "UNRESOLVED",
            "generation_script_or_receipt": "UNRESOLVED",
            "runtime_fetch_revision": "UNPINNED",
            "checkpoint_binding": "UNRESOLVED",
            "biological_admission": "blocked",
        },
        "reference_identity": {
            "publication_tutorial_label": "UCSC_hg38",
            "publication_tutorial_url": UCSC_HG38_FASTA_URL,
            "current_external_checksum_used_for_admission": False,
            "GRCh38_p14_equivalence_claimed": False,
            "exact_checkpoint_reference_identity": "UNRESOLVED",
        },
        "checkpoint": {
            "path": "models/human_model.pt",
            "lfs_sha256": CHECKPOINT_LFS_SHA256,
            "size_bytes": CHECKPOINT_SIZE_BYTES,
            "object_downloaded": False,
            "deserialized": False,
            "state_key_shape_dtype_inventory": "UNRESOLVED",
            "loader_evidence": "model.load_state_dict(torch.load(path,map_location=cpu))",
            "serialization_disposition": "pickle_capable_torch_load_path_not_executed",
            "safe_tensor_conversion": "UNAVAILABLE",
        },
        "terms": {
            "publication_code": "MIT",
            "space_pinned": space_terms,
            "reference_pinned": reference_terms,
            "space_current": api_terms(space_current_api, None),
            "reference_current": api_terms(reference_current_api, None),
            "license_endpoint_http_status": license_status,
            "checkpoint_weights": "UNDECLARED",
            "reference_arrays": "UNDECLARED",
            "derivative_weight_redistribution": "UNRESOLVED",
        },
        "released_execution_path": {
            "entrypoint": "app.py_to_func_gradio.run_epcotv2",
            "source_revision": SPACE_REVISION,
            "source_file_sha256": source_hashes,
            "user_ATAC_input": "pickle.load_UNSAFE",
            "prediction_output": "pickle.dump",
            "checkpoint_load": "torch.load_without_safe_weights_only_UNSAFE",
            "reference_fetch": "hf_hub_download_without_revision_MUTABLE",
            "publication_torch": "1.10.1",
            "space_torch": "2.2.1",
            "space_huggingface_hub_version": "UNPINNED",
            "safe_local_execution_path": "UNAVAILABLE",
            "numeric_parity": "UNTESTED",
        },
        "resolved_items": [
            "publication_tutorial_A_C_G_T_channel_order_and_N_encoding",
            "publication_tutorial_UCSC_hg38_download_route",
            "paper_linked_checkpoint_LFS_identity_and_size",
            "released_Space_entrypoint_and_unsafe_loader_paths",
        ],
        "blocking_items": [
            "Space_reference_NPZ_row_to_base_mapping_not_declared",
            "checkpoint_bound_FASTA_identity_and_checksum_absent",
            "Space_checkpoint_and_reference_array_terms_undeclared",
            "checkpoint_state_inventory_unverified_without_deserialization",
            "released_reference_fetch_is_not_revision_pinned",
            "released_execution_path_accepts_and_emits_pickle",
            "publication_and_Space_runtime_versions_diverge",
            "safe_local_adapter_and_numeric_parity_fixture_absent",
        ],
        "terminal_disposition": {
            "biological_execution": "blocked",
            "checkpoint_download": "blocked_terms",
            "gpu_admission": False,
            "champion_eligibility": "ineligible",
            "required_next_evidence": [
                "written_terms_or_institutional_legal_determination_for_Space_weights_and_reference_arrays",
                "checkpoint_bound_reference_generation_receipt_with_FASTA_checksum_and_A_C_G_T_row_order",
                "safe_checkpoint_conversion_or_verified_weights_only_state_inventory",
                "fully_pinned_local_runtime_and_all_head_numeric_parity_fixture",
            ],
        },
        "checkpoint_object_present_in_audit": False,
        "project_data_read": False,
        "outcomes_read": False,
        "sealed_labels_read": False,
        "upstream_code_imported_or_executed": False,
        "pickle_or_torch_deserialization_performed": False,
    }
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "audit.json").write_text(
        json.dumps(result, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--admission", type=Path, required=True)
    parser.add_argument("--source-cache", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise EPCOTv2ResolutionError(f"output already exists: {args.output}")
    result = audit(args)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
