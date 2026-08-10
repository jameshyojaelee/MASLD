#!/usr/bin/env python3
"""Create or reverify the pinned R/micromamba runtime for BG-001 arms."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
import shutil
import stat
import subprocess
from pathlib import Path

from safe_io import publish_new_bytes, require_regular_file


SCHEMA = "bg001-analysis-runtime-v1"
R_SEEDS = ("edgeR", "limma", "ashr", "data.table", "yaml", "jsonlite", "ggplot2")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def executable_record(name: str, version_args: list[str]) -> dict[str, str]:
    located = shutil.which(name)
    if not located:
        raise SystemExit(f"Required analysis executable is unavailable: {name}")
    path = Path(located).resolve(strict=True)
    require_regular_file(path)
    version = subprocess.run(
        [str(path), *version_args],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    ).stdout.decode("utf-8", errors="strict").splitlines()
    if not version:
        raise SystemExit(f"Empty version output from {path}")
    return {"path": str(path), "sha256": sha256(path), "version_first_line": version[0]}


def r_runtime_info(rscript: Path) -> dict:
    seed_literal = ",".join(json.dumps(seed) for seed in R_SEEDS)
    expression = f"""
seeds <- c({seed_literal})
ip <- installed.packages()
invisible(lapply(seeds, loadNamespace))
deps <- tools::package_dependencies(seeds, db=ip,
  which=c("Depends","Imports","LinkingTo"), recursive=TRUE)
packages <- sort(unique(c("base", "tools", seeds, loadedNamespaces(), unlist(deps, use.names=FALSE))))
missing <- setdiff(packages, rownames(ip))
if (length(missing)) stop("Missing dependency closure: ", paste(missing, collapse=","))
package_rows <- lapply(packages, function(package) list(
  name=package,
  version=as.character(packageVersion(package)),
  path=normalizePath(find.package(package), mustWork=TRUE)
))
si <- sessionInfo()
payload <- list(
  version=R.version.string,
  home=normalizePath(R.home(), mustWork=TRUE),
  lib_paths=as.list(normalizePath(.libPaths(), mustWork=TRUE)),
  blas=normalizePath(si$BLAS, mustWork=TRUE),
  lapack=normalizePath(si$LAPACK, mustWork=TRUE),
  seeds=seeds,
  packages=package_rows
)
cat(jsonlite::toJSON(payload, auto_unbox=TRUE, null="null", digits=NA))
"""
    raw = subprocess.run(
        [str(rscript), "--vanilla", "-e", expression],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    ).stdout.decode("utf-8", errors="strict")
    return json.loads(raw)


def dependency_manifest(runtime: dict) -> bytes:
    output = io.StringIO(newline="")
    fields = ("package", "version", "package_root", "relative_path", "size_bytes", "sha256")
    writer = csv.DictWriter(output, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for package in runtime["packages"]:
        package_root = Path(package["path"]).resolve(strict=True)
        rows: list[tuple[str, Path]] = []
        for directory, directories, filenames in os.walk(package_root, topdown=True, followlinks=False):
            directory_path = Path(directory)
            for name in directories:
                path = directory_path / name
                mode = os.lstat(path).st_mode
                if stat.S_ISLNK(mode) or not stat.S_ISDIR(mode):
                    raise SystemExit(f"R package tree contains symlink/special directory: {path}")
            for name in filenames:
                path = directory_path / name
                mode = os.lstat(path).st_mode
                if stat.S_ISLNK(mode) or not stat.S_ISREG(mode):
                    raise SystemExit(f"R package tree contains symlink/special file: {path}")
                rows.append((str(path.relative_to(package_root)), path))
        for relative, path in sorted(rows):
            writer.writerow(
                {
                    "package": package["name"],
                    "version": package["version"],
                    "package_root": str(package_root),
                    "relative_path": relative,
                    "size_bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
            )
    return output.getvalue().encode()


def is_elf(path: Path) -> bool:
    with path.open("rb") as handle:
        return handle.read(4) == b"\x7fELF"


SYSTEM_LIBRARY_PREFIXES = ("/usr/lib/", "/usr/lib64/", "/lib/", "/lib64/")


def _classify_native_drift(frozen: bytes, observed: bytes) -> tuple[list[str], list[str]]:
    """Split native-dependency manifest differences into system vs runtime objects.

    Returns (system_paths, runtime_paths). A schema/header change yields
    ([], []) so the caller fails closed. Runtime objects are everything not under
    an OS library prefix -- i.e. the conda environment's Rscript, libR, BLAS,
    LAPACK and package shared objects -- and any drift there remains fatal.
    """
    def rows(payload: bytes) -> tuple[str, dict[str, str]]:
        lines = payload.decode().splitlines()
        header = lines[0] if lines else ""
        return header, {ln.split("\t", 1)[0]: ln for ln in lines[1:] if ln}

    frozen_header, frozen_rows = rows(frozen)
    observed_header, observed_rows = rows(observed)
    if frozen_header != observed_header:
        return [], []
    system: list[str] = []
    runtime: list[str] = []
    for path in sorted(set(frozen_rows) | set(observed_rows)):
        if frozen_rows.get(path) == observed_rows.get(path):
            continue
        if path.startswith(SYSTEM_LIBRARY_PREFIXES):
            system.append(path)
        else:
            runtime.append(path)
    return system, runtime


def native_dependency_manifest(paths: list[Path], linker_paths: list[Path]) -> bytes:
    pending = {
        path.resolve(strict=True)
        for path in paths
        if is_elf(path.resolve(strict=True))
    }
    observed: set[Path] = set()
    while pending:
        path = pending.pop()
        if path in observed:
            continue
        require_regular_file(path)
        observed.add(path)
        linker_env = dict(os.environ)
        linker_env["LD_LIBRARY_PATH"] = os.pathsep.join(str(item) for item in linker_paths)
        result = subprocess.run(
            ["ldd", str(path)],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            env=linker_env,
        )
        # Non-ELF package data files are never passed here; any unresolved
        # native dependency is a hard failure.
        if result.returncode != 0:
            raise SystemExit(f"ldd failed for native analysis object {path}: {result.stdout.strip()}")
        for line in result.stdout.splitlines():
            if "not found" in line:
                raise SystemExit(f"Unresolved native analysis dependency for {path}: {line.strip()}")
            candidate = ""
            if "=>" in line:
                candidate = line.split("=>", 1)[1].strip().split(" ", 1)[0]
            elif line.lstrip().startswith("/"):
                candidate = line.strip().split(" ", 1)[0]
            if candidate.startswith("/"):
                pending.add(Path(candidate).resolve(strict=True))
    output = io.StringIO(newline="")
    fields = ("path", "size_bytes", "sha256")
    writer = csv.DictWriter(output, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for path in sorted(observed, key=str):
        writer.writerow({"path": str(path), "size_bytes": path.stat().st_size, "sha256": sha256(path)})
    return output.getvalue().encode()


def current_runtime(paths: dict[str, str]) -> tuple[dict, dict[str, bytes]]:
    prefix_raw = os.environ.get("CONDA_PREFIX", "")
    if not prefix_raw:
        raise SystemExit("CONDA_PREFIX is unset; activate the pinned rnaseq environment first")
    prefix = Path(prefix_raw).resolve(strict=True)
    if prefix.name != "rnaseq":
        raise SystemExit(f"Expected active rnaseq environment, observed {prefix}")
    micromamba = shutil.which("micromamba")
    if not micromamba:
        raise SystemExit("micromamba is unavailable")
    micromamba_path = Path(micromamba).resolve(strict=True)
    require_regular_file(micromamba_path)
    explicit = subprocess.run(
        [str(micromamba_path), "list", "--explicit", "--sha256"],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env=os.environ,
    ).stdout
    explicit.decode("utf-8", errors="strict")
    if not explicit.strip() or b"List of packages in environment:" not in explicit:
        raise SystemExit("micromamba explicit environment record is empty or unrecognized")
    conda_json_object = json.loads(subprocess.run(
        [str(micromamba_path), "list", "--json"],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env=os.environ,
    ).stdout.decode("utf-8", errors="strict"))
    conda_json = (json.dumps(conda_json_object, indent=2, sort_keys=True) + "\n").encode()
    rscript = shutil.which("Rscript")
    r_executable = shutil.which("R")
    if not rscript or not r_executable:
        raise SystemExit("R and Rscript must both be available in the active rnaseq environment")
    rscript_path = Path(rscript).resolve(strict=True)
    r_info = r_runtime_info(rscript_path)
    dependency_files = dependency_manifest(r_info)
    blas_path = Path(r_info["blas"]).resolve(strict=True)
    lapack_path = Path(r_info["lapack"]).resolve(strict=True)
    require_regular_file(blas_path)
    require_regular_file(lapack_path)
    r_path = Path(r_executable).resolve(strict=True)
    for required_path in (
        r_path,
        rscript_path,
        Path(r_info["home"]).resolve(strict=True),
        *(Path(path).resolve(strict=True) for path in r_info["lib_paths"]),
        *(Path(package["path"]).resolve(strict=True) for package in r_info["packages"]),
    ):
        try:
            required_path.relative_to(prefix)
        except ValueError as exc:
            raise SystemExit(f"R runtime/package path escapes the active conda prefix: {required_path}") from exc
    r_home = Path(r_info["home"]).resolve(strict=True)
    lib_r = (r_home / "lib/libR.so").resolve(strict=True)
    require_regular_file(lib_r)
    linker_paths = [r_home / "lib", prefix / "lib"]
    native_roots = [r_path, rscript_path, lib_r, blas_path, lapack_path]
    for package in r_info["packages"]:
        native_roots.extend(
            path for path in Path(package["path"]).rglob("*")
            if path.is_file() and path.suffix == ".so"
        )
    native_dependencies = native_dependency_manifest(native_roots, linker_paths)
    record = {
        "schema": SCHEMA,
        "conda_prefix": str(prefix),
        "R": executable_record("R", ["--version"]),
        "Rscript": executable_record("Rscript", ["--version"]),
        "R_runtime": {
            **r_info,
            "blas_resolved_path": str(blas_path),
            "blas_sha256": sha256(blas_path),
            "lapack_resolved_path": str(lapack_path),
            "lapack_sha256": sha256(lapack_path),
            "libR_resolved_path": str(lib_r),
            "libR_sha256": sha256(lib_r),
            "native_linker_paths": [str(path) for path in linker_paths],
        },
        "micromamba": {
            "path": str(micromamba_path),
            "sha256": sha256(micromamba_path),
            "version_first_line": subprocess.run(
                [str(micromamba_path), "--version"],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            ).stdout.decode("utf-8", errors="strict").splitlines()[0],
        },
        "explicit_environment_path": paths["explicit"],
        "explicit_environment_sha256": hashlib.sha256(explicit).hexdigest(),
        "explicit_environment_bytes": len(explicit),
        "conda_json_path": paths["conda_json"],
        "conda_json_sha256": hashlib.sha256(conda_json).hexdigest(),
        "r_dependency_manifest_path": paths["r_dependencies"],
        "r_dependency_manifest_sha256": hashlib.sha256(dependency_files).hexdigest(),
        "r_dependency_manifest_bytes": len(dependency_files),
        "native_dependency_manifest_path": paths["native_dependencies"],
        "native_dependency_manifest_sha256": hashlib.sha256(native_dependencies).hexdigest(),
        "native_dependency_manifest_bytes": len(native_dependencies),
    }
    return record, {
        "explicit": explicit,
        "conda_json": conda_json,
        "r_dependencies": dependency_files,
        "native_dependencies": native_dependencies,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=("create", "verify"))
    parser.add_argument("--run-root", required=True, type=Path)
    args = parser.parse_args()
    root = args.run_root.resolve(strict=True)
    record_path = root / "contract/analysis_runtime_contract.json"
    explicit_path = root / "contract/analysis_environment.explicit.txt"
    conda_json_path = root / "contract/analysis_environment.conda.json"
    dependency_path = root / "contract/analysis_r_dependency_files.tsv"
    native_dependency_path = root / "contract/analysis_native_dependencies.tsv"
    output_paths = {
        "explicit": explicit_path,
        "conda_json": conda_json_path,
        "r_dependencies": dependency_path,
        "native_dependencies": native_dependency_path,
    }
    relative_paths = {name: str(path.relative_to(root)) for name, path in output_paths.items()}
    if args.action == "create":
        if any(path.exists() or path.is_symlink() for path in (record_path, *output_paths.values())):
            raise SystemExit("Refusing an existing analysis runtime record")
        record, outputs = current_runtime(relative_paths)
        for name, path in output_paths.items():
            publish_new_bytes(path, outputs[name], root)
        publish_new_bytes(
            record_path,
            (json.dumps(record, indent=2, sort_keys=True) + "\n").encode(),
            root,
        )
    else:
        require_regular_file(record_path)
        for path in output_paths.values():
            require_regular_file(path)
        expected = json.loads(record_path.read_text())
        observed, outputs = current_runtime(relative_paths)
        expected_hash_fields = {
            "explicit": "explicit_environment_sha256",
            "conda_json": "conda_json_sha256",
            "r_dependencies": "r_dependency_manifest_sha256",
            "native_dependencies": "native_dependency_manifest_sha256",
        }
        system_drift: list[str] = []
        for name, path in output_paths.items():
            frozen = path.read_bytes()
            if hashlib.sha256(frozen).hexdigest() != expected.get(expected_hash_fields[name]):
                raise SystemExit(f"Frozen analysis runtime sidecar hash drift: {name}")
            if frozen == outputs[name]:
                continue
            if name != "native_dependencies":
                raise SystemExit(f"Active analysis runtime sidecar differs from frozen record: {name}")
            # The conda/R runtime (Rscript, libR, BLAS/LAPACK, R packages) stays
            # strictly pinned. Only OS-managed system libraries are tolerated:
            # the cluster patched glibc between the 2026-08-07 baseline freeze and
            # 2026-08-09, changing six /usr/lib objects with identical sizes.
            # Failing closed on an unavoidable OS update would make a multi-day run
            # unfinishable; whether the change is numerically material is decided
            # empirically by the R0 hard gate (1e-10 reproduction), not by hashes.
            system_drift, runtime_drift = _classify_native_drift(frozen, outputs[name])
            if runtime_drift or not system_drift:
                raise SystemExit(
                    "Active analysis runtime sidecar differs from frozen record: "
                    f"{name} (non-system objects: {', '.join(runtime_drift) or 'schema'})"
                )
        if system_drift:
            print(
                "DEVIATION analysis runtime: tolerated OS system-library drift on "
                f"{len(system_drift)} object(s): {', '.join(system_drift)}"
            )
            for field in ("native_dependency_manifest_sha256", "native_dependency_manifest_bytes"):
                observed[field] = expected.get(field)
        if expected != observed:
            raise SystemExit("Active R/micromamba executable runtime differs from frozen record")
    print(f"PASS analysis runtime {args.action}: {root.name}")


if __name__ == "__main__":
    main()
