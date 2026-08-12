#!/usr/bin/env python3
"""Independently validate the accepted spatial diagnostic bundle index."""

from __future__ import annotations

import argparse
import importlib.util
import re
import subprocess
from pathlib import Path

import pandas as pd


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    "spatial_diagnostic_bundle_builder", HERE / "45_build_spatial_diagnostic_bundle.py"
)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("cannot import spatial diagnostic bundle builder")
builder = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(builder)


def pdf_pages(path: Path) -> int:
    text = subprocess.run(["pdfinfo", str(path)], check=True, capture_output=True, text=True).stdout
    for line in text.splitlines():
        if line.startswith("Pages:"):
            return int(line.split(":", 1)[1].strip())
    raise RuntimeError(f"pdfinfo did not report pages: {path}")


def validate(project: Path, output: Path) -> None:
    project = project.resolve()
    output = output.resolve()
    if not (output / "BUILD_COMPLETE").is_file() or (output / "READY").exists():
        raise RuntimeError("bundle is absent, incomplete, or already validated")
    registry = pd.read_csv(output / "diagnostic_registry.tsv", sep="\t", dtype=str)
    artifact = pd.read_csv(output / "artifact_manifest.tsv", sep="\t", dtype=str)
    code = pd.read_csv(output / "code_manifest.tsv", sep="\t", dtype=str)
    expected_ids = {component["component_id"] for component in builder.COMPONENTS}
    if len(registry) != 5 or set(registry["component_id"]) != expected_ids:
        raise RuntimeError("diagnostic registry is not the complete five-component family")
    if set(registry["canonical_write_allowed"]) != {"FALSE"}:
        raise RuntimeError("diagnostic bundle authorizes canonical writes")
    if int(registry["n_pdfs"].astype(int).sum()) != 6:
        raise RuntimeError("diagnostic bundle must contain six PDFs")
    component_map = {component["component_id"]: component for component in builder.COMPONENTS}
    expected_artifacts = []
    for row in registry.itertuples(index=False):
        component = component_map[row.component_id]
        root = project / row.relative_root
        ready = root / "READY"
        if (
            row.candidate_release_id != component["candidate"]
            or builder.read_ready_status(ready) != component["status"]
            or builder.sha256_file(ready) != row.ready_sha256
        ):
            raise RuntimeError(f"component READY drift: {row.component_id}")
        if any(token in row.candidate_release_id for token in ("failed", "superseded")):
            raise RuntimeError("failed or superseded candidate entered the bundle")
        for path in sorted(item for item in root.rglob("*") if item.is_file()):
            expected_artifacts.append((row.component_id, path.relative_to(project).as_posix()))
    observed_artifacts = list(zip(artifact["component_id"], artifact["relative_path"], strict=True))
    if observed_artifacts != expected_artifacts or artifact["relative_path"].duplicated().any():
        raise RuntimeError("artifact manifest path family drift")
    for row in artifact.itertuples(index=False):
        path = project / row.relative_path
        if not path.is_file() or path.stat().st_size != int(row.bytes) or builder.sha256_file(path) != row.sha256:
            raise RuntimeError(f"diagnostic artifact drift: {row.relative_path}")
        if "/figures/main/" in f"/{row.relative_path}" or "/results/" in f"/{row.relative_path}":
            raise RuntimeError(f"canonical artifact leaked into bundle: {row.relative_path}")
    expected_code = builder.code_paths(HERE, project)
    if code["relative_path"].tolist() != [path.relative_to(project).as_posix() for path in expected_code]:
        raise RuntimeError("code manifest path family drift")
    for row in code.itertuples(index=False):
        path = project / row.relative_path
        if not path.is_file() or path.stat().st_size != int(row.bytes) or builder.sha256_file(path) != row.sha256:
            raise RuntimeError(f"diagnostic code drift: {row.relative_path}")
    pdf_rows = artifact[artifact["relative_path"].str.endswith(".pdf")]
    if len(pdf_rows) != 6:
        raise RuntimeError("artifact manifest PDF count drift")
    for row in pdf_rows.itertuples(index=False):
        path = project / row.relative_path
        if pdf_pages(path) != 1 or re.search(rb"/Subtype\s*/Type3\b", path.read_bytes()):
            raise RuntimeError(f"invalid diagnostic PDF: {row.relative_path}")
    build = pd.read_csv(output / "BUILD_COMPLETE", sep="\t", dtype=str)
    if len(build) != 1 or build.iloc[0]["status"] != "built_pending_independent_validation":
        raise RuntimeError("invalid diagnostic bundle build marker")
    expected_hashes = {
        "diagnostic_registry_sha256": builder.sha256_file(output / "diagnostic_registry.tsv"),
        "artifact_manifest_sha256": builder.sha256_file(output / "artifact_manifest.tsv"),
        "code_manifest_sha256": builder.sha256_file(output / "code_manifest.tsv"),
    }
    if any(build.iloc[0][key] != value for key, value in expected_hashes.items()):
        raise RuntimeError("bundle build hashes drifted")
    pd.DataFrame([{
        "bundle_release_id": builder.RELEASE_ID,
        "status": "validated_spatial_diagnostic_bundle_awaiting_adjudication",
        "n_components": 5,
        "n_pdfs": 6,
        "n_artifacts": len(artifact),
        "n_code_files": len(code),
        **expected_hashes,
        "canonical_write_allowed": "FALSE",
    }]).to_csv(output / "READY", sep="\t", index=False)
    print(f"PASS spatial diagnostic bundle: 5 components, 6 PDFs, {len(artifact)} artifacts")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    validate(args.project_root, args.output_root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
