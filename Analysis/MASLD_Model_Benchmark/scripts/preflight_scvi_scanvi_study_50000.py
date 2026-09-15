#!/usr/bin/env python3
"""Fail-closed inclusion audit for the 50k scVI/scANVI campaign."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
from importlib import metadata as importlib_metadata
import json
from pathlib import Path
from typing import Any

from masld_bench.artifacts import verify_frozen_tree


EXPECTED_SOURCE_SHA256 = "672b1f75a1c3ca264066dd2a0bf24e479397d8ec5e2f67f13ccdd4be14879509"
EXPECTED_SPLIT_SHA256 = "10927e162581375d866adbcf2b3bb7bd0dc4fda9fd6dfe033198df7d95be2680"
EXPECTED_WHEEL_SHA256 = "c453b75b0fa0d222bf0eb3a531607d17f5ea0dbe5b848ab12fd8cc24599884e5"
EXPECTED_CODE_COMMIT = "dc95497"
EXPECTED_SCVI_LICENSE_SHA256 = "66399db0284d2539790efb348886ab0c1f745bbe5fea6ac38a00465a14adc8f5"
EXPECTED_SEEDS = (20260824, 20260825, 20260826)


class SCVIPreflightError(ValueError):
    """Raised when the read-only campaign is not production-eligible."""


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise SCVIPreflightError(f"TSV header missing: {path}")
        return list(reader)


def run(
    *,
    source: Path,
    split: Path,
    runtime: Path,
    wrapper: Path,
    output: Path,
) -> None:
    if output.exists():
        raise SCVIPreflightError("refusing to overwrite preflight output")
    if sha256_file(source / "ARTIFACTS.json") != EXPECTED_SOURCE_SHA256:
        raise SCVIPreflightError("source artifact hash differs")
    if sha256_file(split / "ARTIFACTS.json") != EXPECTED_SPLIT_SHA256:
        raise SCVIPreflightError("split artifact hash differs")
    source_artifacts = verify_frozen_tree(source)
    split_artifacts = verify_frozen_tree(split)
    runtime_artifacts = verify_frozen_tree(runtime)
    source_manifest = json.loads(
        (source / "frozen_screen_manifest.json").read_text(encoding="utf-8")
    )
    runtime_metadata = runtime_artifacts["metadata"]
    if (
        source_artifacts["metadata"].get("row_count") != 50_000
        or source_manifest.get("input_contract", {}).get("matrix")
        != "X_raw_nonnegative_integer_UMI_counts"
        or source_manifest.get("matrix", {}).get("raw_umi_counts_preserved") is not True
        or source_manifest.get("matrix", {}).get("shape") != [50_000, 37_533]
        or source_manifest.get("matrix", {}).get("nnz") != 76_363_854
        or source_manifest.get("selection", {}).get("sealed_outcomes_used") is not False
        or split_artifacts["metadata"].get("donors") != 102
        or split_artifacts["metadata"].get("studies") != 7
        or split_artifacts["metadata"].get("target_labels_used_for_assignment") is not False
        or runtime_metadata.get("scvi_tools_version") != "1.4.2"
        or runtime_metadata.get("scvi_tools_wheel_sha256") != EXPECTED_WHEEL_SHA256
        or runtime_metadata.get("code_commit") != EXPECTED_CODE_COMMIT
        or runtime_metadata.get("portable_site_packages") is not True
        or runtime_metadata.get("direct_inductive_query_probe") != "passed"
        or runtime_metadata.get("query_training_allowed") is not False
        or runtime_metadata.get("query_adaptation_allowed") is not False
        or runtime_metadata.get("project_data_read") is not False
        or runtime_metadata.get("sealed_outcomes_read") is not False
    ):
        raise SCVIPreflightError("dataset, split, or runtime admission differs")

    split_rows = _rows(split / "row_outer_folds.tsv")
    inner_rows = _rows(split / "inner_donor_folds.tsv")
    if len(split_rows) != 50_000:
        raise SCVIPreflightError("outer split row count differs")
    donor_outer: dict[str, set[int]] = {}
    study_outer: dict[str, set[int]] = {}
    for row in split_rows:
        outer = int(row["outer_fold"])
        donor_outer.setdefault(row["donor_id"], set()).add(outer)
        study_outer.setdefault(row["dataset"], set()).add(outer)
    if (
        set().union(*donor_outer.values()) != set(range(5))
        or any(len(values) != 1 for values in donor_outer.values())
        or any(len(values) != 1 for values in study_outer.values())
    ):
        raise SCVIPreflightError("donor or study crosses outer folds")
    inner_pairs = {(int(row["held_outer_fold"]), row["donor_id"]) for row in inner_rows}
    if len(inner_pairs) != len(inner_rows):
        raise SCVIPreflightError("inner donor assignment repeats")
    for held in range(5):
        training_donors = {
            donor for donor, folds in donor_outer.items() if held not in folds
        }
        observed = {
            row["donor_id"] for row in inner_rows if int(row["held_outer_fold"]) == held
        }
        observed_folds = {
            int(row["inner_fold"])
            for row in inner_rows
            if int(row["held_outer_fold"]) == held
        }
        if observed != training_donors or observed_folds != set(range(5)):
            raise SCVIPreflightError("inner donor roster differs")

    wrapper_text = wrapper.read_text(encoding="utf-8")
    script_text = (
        wrapper.parents[1] / "scripts/fit_predict_scvi_scanvi_study_50000.py"
    ).read_text(encoding="utf-8")
    header = "\n".join(
        line for line in wrapper_text.splitlines() if line.startswith("#SBATCH")
    )
    if (
        "--partition=gpu" not in header
        or "--account=nslab" not in header
        or "--qos=nslab" not in header
        or "--gres=gpu:l40s:1" not in header.lower()
        or "--array" in header
        or "load_query_data" in script_text
        or "prepare_query_anndata" in script_text
        or "LogisticRegression" in script_text
        or "labels[query]" in script_text
        or "early_stopping=False" not in script_text
        or "train_size=1.0" not in script_text
        or "MASLD_REQUIRE_CUDA=1" not in wrapper_text
    ):
        raise SCVIPreflightError("GPU wrapper or query firewall differs")

    output.mkdir(mode=0o750)
    dependency_licenses: list[dict[str, Any]] = []
    for distribution in sorted(
        importlib_metadata.distributions(),
        key=lambda value: (value.metadata.get("Name", "").lower(), value.version),
    ):
        package_metadata = distribution.metadata
        dependency_licenses.append(
            {
                "name": package_metadata.get("Name", "UNRESOLVED"),
                "version": distribution.version,
                "license": package_metadata.get("License", "UNRESOLVED"),
                "license_expression": package_metadata.get(
                    "License-Expression", "UNRESOLVED"
                ),
                "license_classifiers": sorted(
                    value
                    for value in package_metadata.get_all("Classifier", [])
                    if value.startswith("License ::")
                ),
                "license_files": sorted(
                    package_metadata.get_all("License-File", [])
                ),
            }
        )
    scvi_distribution = importlib_metadata.distribution("scvi-tools")
    scvi_license = scvi_distribution.locate_file(
        "scvi_tools-1.4.2.dist-info/licenses/LICENSE"
    )
    if (
        scvi_distribution.version != "1.4.2"
        or not scvi_distribution.metadata.get("License", "").startswith(
            "BSD 3-Clause License"
        )
        or sha256_file(Path(scvi_license)) != EXPECTED_SCVI_LICENSE_SHA256
    ):
        raise SCVIPreflightError("scvi-tools code license authority differs")
    with (output / "dependency_license_census.json").open(
        "x", encoding="utf-8"
    ) as handle:
        json.dump(
            {
                "schema_version": "masld-bench-python-license-census-v1",
                "scvi_tools_license": "BSD-3-Clause",
                "scvi_tools_license_sha256": EXPECTED_SCVI_LICENSE_SHA256,
                "dependencies": dependency_licenses,
                "unresolved_dependency_terms_require_release_review": any(
                    record["license"] == "UNRESOLVED"
                    and record["license_expression"] == "UNRESOLVED"
                    and not record["license_classifiers"]
                    and not record["license_files"]
                    for record in dependency_licenses
                ),
            },
            handle,
            sort_keys=True,
            separators=(",", ":"),
        )
        handle.write("\n")
    payload: dict[str, Any] = {
        "schema_version": "masld-bench-scvi-scanvi-preflight-v1",
        "status": "pass_build_only_gpu_pending",
        "dataset_view_id": "resource_atlas_frozen_screen_50000_v1",
        "split_id": "resource_atlas_study_outer_5fold_v1",
        "rows": 50_000,
        "donors": 102,
        "studies": 7,
        "outer_split_unit": "study",
        "inner_split_unit": "donor",
        "seeds": list(EXPECTED_SEEDS),
        "runtime_artifacts_sha256": sha256_file(runtime / "ARTIFACTS.json"),
        "runtime_direct_query_probe": "passed_on_cpu",
        "scvi_tools_code_license": "BSD-3-Clause",
        "dependency_license_census_written": True,
        "dependency_terms_release_review_required": True,
        "gpu_probe_pending": True,
        "query_training_allowed": False,
        "query_adaptation_allowed": False,
        "common_head_fit": False,
        "project_data_model_fit_executed": False,
        "metrics_calculated": False,
        "sealed_outcomes_read": False,
    }
    with (output / "preflight.json").open("x", encoding="utf-8") as handle:
        json.dump(payload, handle, sort_keys=True, separators=(",", ":"))
        handle.write("\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--split", required=True, type=Path)
    parser.add_argument("--runtime", required=True, type=Path)
    parser.add_argument("--wrapper", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    run(
        source=arguments.source,
        split=arguments.split,
        runtime=arguments.runtime,
        wrapper=arguments.wrapper,
        output=arguments.output,
    )
    print(json.dumps({"output": str(arguments.output.resolve()), "status": "passed"}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
