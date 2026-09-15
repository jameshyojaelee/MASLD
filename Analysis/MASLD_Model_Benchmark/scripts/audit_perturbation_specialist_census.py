#!/usr/bin/env python3
"""Audit perturbation specialists, native task boundaries, and replication fail closed."""

from __future__ import annotations

import argparse
import hashlib
import json
import tomllib
from pathlib import Path
from typing import Any

from masld_bench.artifacts import ArtifactError, verify_frozen_tree


class CensusAuditError(RuntimeError):
    """Raised when a frozen perturbation authority or requirement differs."""


CENSUS_RELATIVE = (
    "config/artifacts/model-authorities/perturbation_specialists/census.json"
)
CENSUS_SHA256 = "5cd29369a6b3691bd428ab218fccec9c3f800420c8bcaa994603b663aa2c2450"

AUTHORITY_SHA256 = {
    "config/artifacts/models/gears/checkpoints.json": "221a79988ae00d76d9b6c84ffd85082af5cc3ca710d904883b86e6a695837285",
    "config/artifacts/models/gears/development_crosswalk.json": "366d4277ffbf3e40f058bbae757d9d84d374f35c495bd2be320faa66a8ed4451",
    "config/artifacts/models/gears/exposure_audit.json": "807a2d150875879ca359c458ebcaf723b18a3e264cc80024147ca08017bfc955",
    "config/artifacts/models/scgen/checkpoints.json": "95347a1f30eaddee5286c6ef865290807d5ed940a2240b966b04c235f2acf124",
    "config/artifacts/models/scgen/development_crosswalk.json": "352ae2e6b5603b1f4d16510b5de8dfe52dde47af0f1f231d486624acf46c29c8",
    "config/artifacts/models/scgen/exposure_audit.json": "06fd2d03467a73861426cc12134ec2596e756b778cb18087898ad478066cdce7",
    "config/artifacts/models/cpa/checkpoints.json": "99af16654c9042b34f564e7ee8c797be66d658f163aeb19744ac5ac448ab7c1f",
    "config/artifacts/models/cpa/development_crosswalk.json": "41777bf4e42ae19adb2bfb5d75ae0cc1be7a3ee9066263f5a139b23f5a30d265",
    "config/artifacts/models/cpa/exposure_audit.json": "625005ad12c9c2edf3a35f9c80dad5aab7ff1973aad9f3eb9247819851a06a32",
    "config/artifacts/models/cellot/checkpoints.json": "c6ab79d62c4595eaad2ca81580427e1c9fb05bc51bdcd61c618c7bdd4c4e9261",
    "config/artifacts/models/cellot/development_crosswalk.json": "3e23b10ef1912e69701b1d01ec7991c48ab14b14797158625cf84131c25a851a",
    "config/artifacts/models/cellot/exposure_audit.json": "baed1e149ff3fdbe8ccfe6f2bdd35cd7c3941ce2f06a140660d559736c8187e2",
    "config/artifacts/models/genepert/checkpoints.json": "936774371210db4f5c19f327dbbf3b3185431fa14542fbb242fe3153717fa062",
    "config/artifacts/models/genepert/development_crosswalk.json": "335c930c9fe5defa4801f809fd5a62ca05bd79cc60d2e3edcf8035c265c7d565",
    "config/artifacts/models/genepert/exposure_audit.json": "3bb3761e7281be11ab88afef0709ee1908526934443a8e2632f39fc18251c88c",
    "config/artifacts/models/scgpt_perturbation_head/checkpoints.json": "0d6b1b7c981bf5d180eb6d715d1a2e5b33f4c09464c26ade309e8c9004b01023",
    "config/artifacts/models/scgpt_perturbation_head/development_crosswalk.json": "eeabefa4c0894462d0bfb1fe8acba34c8e1297c65d2417736f797572ce26cf17",
    "config/artifacts/models/scgpt_perturbation_head/exposure_audit.json": "e9735d361801ee4e0d16bb1a70fdb8e535ea7d941180ff301f5420b74c10e859",
    "config/artifacts/models/scfoundation/checkpoints.json": "e56ad920d44ed81b5903df4843925b0442b69df358acaa1b12892de860c0cfc1",
    "config/artifacts/models/scfoundation/development_crosswalk.json": "5be876f9cffb14588af89891370ed09cb985fa7a81ec9da05bc719feed2673aa",
    "config/artifacts/models/scfoundation/exposure_audit.json": "552fd7f8d12ed98d6619f411926ce938290a75309f856689f9364724ed2559cc",
}

CONFIG_SHA256 = {
    "config/models/perturbation.toml": "30afc3d28d6ecae4f3924eee0a149e0a55859514ee9c5dbb54c79c847109e10d",
    "config/models/restricted_unclear.toml": "1bbe758b19e13989c795ec1d1c444523ee4a50e89ea99823194d5458fad4c364",
    "config/tasks/perturbation_transfer.toml": "c44ab6f4b19c54dad8df6a498aa5a2f34bc40811eda6afd7e8bd6745e5574a05",
    "config/tasks/gene_perturbation_response.toml": "a8a0345ae2409df68cf7b30946fae3c45b7e4bd4bd73bfad4b3008e68cac0912",
    "config/splits/gene_perturbation_outer.toml": "279c5be229a61e028914fea741eb74bbd8fd8bf70213728ee82519668e7d8db7",
    "config/datasets/gse281160.toml": "34f6f960a4a734da463a4099359cbc308268d13d7ef49ac89084d3522874d614",
    "config/datasets/gse238219.toml": "96b70b5f781daa3434422d7c936af04c893739f6016847ff70f510f87a2e7f02",
    "config/datasets/gse264667.toml": "178cc128fb9bf267bfc4bef85fab152e310c7ee2268b5ffff4c4f8c33d49d6cf",
    "config/datasets/gse242934.toml": "80d93da8947d03653abbd41aade0f990e7089e2212d1cba8f61a3208f15e443d",
    "config/datasets/cra009621_oscar.toml": "83268a0c950e0c4ca85bcfc6d628db17c3ea9fa3b63d32d20792213600c45e65",
    "config/datasets/gse275483.toml": "12638c09a120228520335abde6bd59e7e513bd04f31f5f25e02860a1f8b5d43a",
    "config/datasets/gse313774.toml": "d90f6a133abb9707df18b04871352d67eeb745c78e987c36a5aa9932cfc919cf",
}

ARTIFACT_SHA256 = {
    "executions/gse238219-target-splits-21067592": "6d6bb5ad785b5aa65573db5e1e1eb15f430bd98cfd163c0571b1fb895a6b31c9",
    "executions/gse264667-h5ad-schema-21067942": "8f6b841073cbd1b113f0257ad3936b8465c97fde349881e3c6e6500a23021fc1",
    "executions/gse264667-metadata-topology-21068136": "2046e8881daa9b7716bd73fb2fdb5d631e273cef4342801c9b2ab872d18e4fd4",
    "executions/gse264667-matrix-scale-21068319": "59708920ad56b388c06722761a4f06a3544b48404fddea47b2309a3026f59723",
    "executions/gse242934-pseudobulk-21067929": "018114c8922cf6b098a1d00ad374312f01cf9f13b5eaeea2ba56162f606d2737",
    "executions/oscar-source-topology-21068446": "809c2cda878f4e53b214a0994d197dab43f981cf4e3a809267919c9839c714ee",
    "executions/perturbation-runtime-gears-21066692_0": "66ac997ddd2ab69ad9a9a23142a8a7488950a10a8bc86315da370a273362c7d4",
    "executions/perturbation-runtime-scgen-21066692_2": "a188a93e6a25a39852965c2a8883a2e17517cc0160d805c3fb0465e56e5086b5",
    "executions/perturbation-runtime-cpa-21066692_3": "9f1100b8049fb5747a54bad42e35c625882c370299a6b83c935ad2644e7f724c",
    "executions/perturbation-runtime-cellot-21066692_4": "f7ef6a30251ed372eb6dbb918c0205fd645f9066bc80f524babc9e425fb65e25",
    "executions/perturbation-runtime-genepert-21066692_5": "492e6af7af5dab623ef0d5a3fd4f6fb6d02a5338d51ebf01632ff6ef2f71993b",
    "executions/perturbation-runtime-scgpt-21066692_1": "5f1f7e719f08714fca1e4e7242d3d4786c6f0cb7bd40881a1abce1ea3de1f50b",
    "executions/gears-core-forward-21067574": "d8ea9af228c89576cfa77361fc1cdc5643d16ec8ab4bd6294da06d26df947c7a",
}


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise CensusAuditError(f"JSON authority is unreadable: {path}: {error}") from error
    if not isinstance(value, dict):
        raise CensusAuditError(f"expected JSON object: {path}")
    return value


def load_toml(path: Path) -> dict[str, Any]:
    try:
        value = tomllib.loads(path.read_text(encoding="utf-8"))
    except (OSError, tomllib.TOMLDecodeError) as error:
        raise CensusAuditError(f"TOML authority is unreadable: {path}: {error}") from error
    if not isinstance(value, dict):
        raise CensusAuditError(f"expected TOML table: {path}")
    return value


def require(condition: bool, message: str) -> None:
    if not condition:
        raise CensusAuditError(message)


def require_hash(root: Path, relative: str, expected: str) -> Path:
    path = root / relative
    if not path.is_file() or digest(path) != expected:
        raise CensusAuditError(f"authority identity differs: {relative}")
    return path


def require_artifact(root: Path, relative: str, expected: str) -> dict[str, Any]:
    path = root / relative
    try:
        verified = verify_frozen_tree(path)
    except ArtifactError as error:
        raise CensusAuditError(f"frozen artifact differs: {relative}: {error}") from error
    if digest(path / "ARTIFACTS.json") != expected:
        raise CensusAuditError(f"artifact manifest identity differs: {relative}")
    return dict(verified)


def model_rows(config: dict[str, Any]) -> dict[str, dict[str, Any]]:
    rows = config.get("models")
    if not isinstance(rows, list):
        raise CensusAuditError("model registry rows differ")
    indexed = {
        str(row.get("model_id")): row for row in rows if isinstance(row, dict)
    }
    if len(indexed) != len(rows):
        raise CensusAuditError("model registry contains duplicate or malformed IDs")
    return indexed


def validate_census(census: dict[str, Any]) -> tuple[set[str], dict[str, dict[str, Any]]]:
    """Validate the central scientific requirements independently of filesystem evidence."""
    expected_models = {
        "gears",
        "scgen",
        "cpa",
        "cellot",
        "genepert",
        "scgpt_perturbation_head",
        "scfoundation_heads",
    }
    census_models = census.get("models")
    require(isinstance(census_models, list), "census model rows differ")
    by_model = {
        str(row.get("model_id")): row
        for row in census_models
        if isinstance(row, dict)
    }
    require(set(by_model) == expected_models, "census model membership differs")
    require(len(by_model) == len(census_models), "duplicate census model IDs")
    for model_id, row in by_model.items():
        require(row.get("runtime_admitted") is False, f"{model_id} runtime was silently admitted")
        require(row.get("schedulable") is False, f"{model_id} became schedulable without a new census")
        require(row.get("open_champion_eligible") is False, f"{model_id} became champion eligible")

    scope = census.get("scope", {})
    require(scope.get("sealed_sources") == [], "perturbation census acquired a sealed source")
    require(scope.get("champion_eligible") is False, "perturbation census became champion eligible")
    require(scope.get("universal_perturbation_objective") is False,
            "native perturbation endpoints were collapsed into one objective")
    global_contract = census.get("global_contract", {})
    forbidden_replicates = set(global_contract.get("not_independent_replicates", []))
    require({"cells", "guides", "target_genes", "sequencing_lanes", "GEM_groups"}
            <= forbidden_replicates, "pseudoreplication guard differs")
    require(global_contract.get("production_jobs_authorized_by_this_census") is False,
            "census must not authorize production jobs")
    require(global_contract.get("sealed_outcomes_accessed") is False,
            "census reports sealed outcome access")
    return expected_models, by_model


def audit(root: Path) -> dict[str, Any]:
    census = load_json(require_hash(root, CENSUS_RELATIVE, CENSUS_SHA256))
    authorities = {
        relative: load_json(require_hash(root, relative, expected))
        for relative, expected in AUTHORITY_SHA256.items()
    }
    configs = {
        relative: load_toml(require_hash(root, relative, expected))
        for relative, expected in CONFIG_SHA256.items()
    }
    artifacts = {
        relative: require_artifact(root, relative, expected)
        for relative, expected in ARTIFACT_SHA256.items()
    }

    expected_models, by_model = validate_census(census)

    perturb_registry = model_rows(configs["config/models/perturbation.toml"])
    for model_id in expected_models - {"scfoundation_heads"}:
        row = perturb_registry.get(model_id, {})
        require(row.get("status") == "deferred" and row.get("admission_blocking") is True,
                f"{model_id} registry is no longer fail closed")
        require(row.get("supported_tasks") == ["gene_perturbation_response"],
                f"{model_id} native task registration differs")
    restricted_registry = model_rows(configs["config/models/restricted_unclear.toml"])
    scfoundation = restricted_registry.get("scfoundation_heads", {})
    require(scfoundation.get("admission_blocking") is True,
            "scFoundation restriction was removed")
    require("gene_perturbation_response" not in scfoundation.get("supported_tasks", []),
            "scFoundation perturbation head was admitted without a separate manifest")

    repository_expectations = {
        "gears": "snap-stanford/GEARS",
        "scgen": "theislab/scgen",
        "cpa": "theislab/cpa",
        "cellot": "bunnech/cellot",
        "genepert": "zou-group/GenePert",
        "scgpt_perturbation_head": "bowang-lab/scGPT",
        "scfoundation_heads": "biomap-research/scFoundation",
    }
    bundle_names = {
        **{model_id: model_id for model_id in expected_models},
        "scfoundation_heads": "scfoundation",
    }
    for model_id, repository in repository_expectations.items():
        checkpoint = authorities[
            f"config/artifacts/models/{bundle_names[model_id]}/checkpoints.json"
        ]
        require(checkpoint.get("repository") == repository,
                f"{model_id} repository identity differs")
    require("UNRESOLVED" in authorities[
        "config/artifacts/models/gears/checkpoints.json"
    ]["terms_audit"]["external_GO_pickle_terms"], "GEARS external prior was silently admitted")
    require(authorities[
        "config/artifacts/models/scfoundation/checkpoints.json"
    ]["task_fit"]["perturbation_response"].startswith("not_admitted"),
            "scFoundation perturbation task disposition differs")

    gene_task = configs["config/tasks/gene_perturbation_response.toml"]
    transfer_task = configs["config/tasks/perturbation_transfer.toml"]
    require(gene_task.get("unit_of_inference") == "independent_culture_replicate_or_animal",
            "gene perturbation inferential unit differs")
    require(gene_task.get("datasets_sealed") == [] and transfer_task.get("datasets_sealed") == [],
            "perturbation task acquired an unreviewed sealed source")
    require(gene_task.get("primary_evaluator_id") == "not_applicable_nonpromotable_v1"
            and transfer_task.get("primary_evaluator_id") == "not_applicable_nonpromotable_v1",
            "a nonpromotable task acquired a champion evaluator")
    require("gse281160" not in gene_task.get("datasets_train", [])
            and "gse281160" not in gene_task.get("datasets_development", []),
            "noncoding GSE281160 entered the gene-target task")

    datasets = census.get("datasets")
    require(isinstance(datasets, list), "census dataset rows differ")
    by_dataset = {
        str(row.get("dataset_id")): row
        for row in datasets
        if isinstance(row, dict)
    }
    require(set(by_dataset) == {
        "gse281160", "gse238219", "gse264667", "gse242934",
        "cra009621_oscar", "gse275483", "gse313774",
    }, "census dataset membership differs")
    require(all(row.get("champion_eligible") is False for row in by_dataset.values()),
            "a current perturbation dataset became champion eligible")
    require(by_dataset["gse264667"]["biological_replication"] == "one_transduced_pool",
            "GSE264667 replication differs")
    require(by_dataset["gse242934"]["biological_replication"] == "one_primary_hepatocyte_pool",
            "GSE242934 replication differs")
    require(by_dataset["gse275483"]["biological_replication"].startswith("four_animals"),
            "GSE275483 replication differs")
    require(by_dataset["gse313774"]["biological_replication"]
            == "three_independent_liver_chip_batches", "GSE313774 replication differs")

    target_splits = load_json(
        root / "executions/gse238219-target-splits-21067592/target_split_manifest.json"
    )
    require(target_splits.get("fold_count") == 11
            and target_splits.get("biological_unit_status") == "UNRESOLVED"
            and target_splits.get("expression_values_used_to_assign_roles") is False,
            "GSE238219 split receipt differs")
    runtime_roots = {
        "gears": "executions/perturbation-runtime-gears-21066692_0",
        "scgen": "executions/perturbation-runtime-scgen-21066692_2",
        "cpa": "executions/perturbation-runtime-cpa-21066692_3",
        "cellot": "executions/perturbation-runtime-cellot-21066692_4",
        "genepert": "executions/perturbation-runtime-genepert-21066692_5",
        "scgpt": "executions/perturbation-runtime-scgpt-21066692_1",
    }
    runtime_status = {}
    for lane, relative in runtime_roots.items():
        receipt = load_json(root / relative / "runtime_audit.json")
        require(receipt.get("status") == "observed_not_admitted"
                and receipt.get("forward_fixture_performed") is False,
                f"{lane} mutable runtime evidence was overclaimed")
        runtime_status[lane] = {
            "primary_imported": receipt.get("primary_imported"),
            "forward_fixture_performed": False,
            "status": "observed_not_admitted",
        }
    require(runtime_status["gears"]["primary_imported"] is True,
            "GEARS import census differs")
    require(all(runtime_status[lane]["primary_imported"] is False
                for lane in {"scgen", "cpa", "cellot", "genepert", "scgpt"}),
            "failed perturbation runtime census differs")
    gears_core = load_json(
        root / "executions/gears-core-forward-21067574/gears_core_forward.json"
    )
    require(gears_core.get("deterministic_exact_repeat") is True
            and gears_core.get("all_finite") is True
            and gears_core.get("biological_data_used") is False
            and gears_core.get("scientific_admission")
            == "blocked; core tensor fixture only",
            "GEARS synthetic core evidence was overclaimed")

    return {
        "schema_version": "masld-bench-perturbation-specialist-census-audit-v1",
        "status": "pass_fail_closed",
        "census_sha256": CENSUS_SHA256,
        "authority_sha256": AUTHORITY_SHA256,
        "config_sha256": CONFIG_SHA256,
        "artifact_manifest_sha256": ARTIFACT_SHA256,
        "models_reconciled": sorted(expected_models),
        "datasets_reconciled": sorted(by_dataset),
        "native_subendpoints_preserved": True,
        "biological_replication_preserved": True,
        "sealed_sources": [],
        "champion_eligible": False,
        "production_jobs_authorized": False,
        "runtime_status": runtime_status,
        "artifact_metadata_verified": sorted(artifacts),
        "terminal_disposition": "no_perturbation_specialist_currently_schedulable",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.project_root.resolve(strict=True)
    if args.output.exists():
        raise CensusAuditError("output already exists")
    result = audit(root)
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "census_audit.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
