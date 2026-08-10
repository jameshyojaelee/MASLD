#!/usr/bin/env python3
"""Capture existing Python/R environments and execution identity without installing anything."""

from __future__ import annotations

import datetime as dt
import hashlib
import importlib.metadata
import os
import platform
import subprocess
import sys

from public_functional_common import CANDIDATE_ROOT, PROJECT_ROOT, sha256_file, write_tsv


def version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not_installed"


def main() -> None:
    rows = [
        {"environment": "rapids_singlecell", "component": "python", "version": platform.python_version()},
        {"environment": "rapids_singlecell", "component": "platform", "version": platform.platform()},
    ]
    for package in ["pandas", "numpy", "scipy", "anndata", "h5py", "openpyxl"]:
        rows.append({"environment": "rapids_singlecell", "component": package, "version": version(package)})
    r_expression = "cat(R.version.string,'\\n');for(p in c('limma','edgeR','readxl','AnnotationDbi','org.Hs.eg.db','RSQLite'))cat(p,as.character(packageVersion(p)),'\\n')"
    output = subprocess.check_output(["micromamba", "run", "-n", "rnaseq", "Rscript", "-e", r_expression], text=True)
    for index, line in enumerate(output.strip().splitlines()):
        fields = line.split(maxsplit=1)
        rows.append({"environment": "rnaseq", "component": "R" if index == 0 else fields[0], "version": fields[-1]})
    write_tsv(CANDIDATE_ROOT / "environment_manifest.tsv", rows, ["environment", "component", "version"])
    git_diff = subprocess.check_output(["git", "diff", "--binary"], cwd=PROJECT_ROOT)
    execution = [{
        "captured_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID", ""),
        "slurm_partition": os.environ.get("SLURM_JOB_PARTITION", ""),
        "git_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=PROJECT_ROOT, text=True).strip(),
        "git_diff_sha256": hashlib.sha256(git_diff).hexdigest(),
        "method_compliance_code_ready_sha256": sha256_file(CANDIDATE_ROOT / "METHOD_COMPLIANCE_EXECUTION_READY.json"),
    }]
    write_tsv(CANDIDATE_ROOT / "execution_manifest.tsv", execution, list(execution[0]))


if __name__ == "__main__":
    main()
