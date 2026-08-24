#!/usr/bin/env python3
"""Acquire bounded primary metadata for liver-relevant gene perturbation sources."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path, PurePosixPath
import tarfile
from typing import Any
from urllib.error import HTTPError
from urllib.request import Request, urlopen


SOURCES: dict[str, dict[str, Any]] = {
    "gse238219": {
        "accessions": ["GSE238219", "PRJNA998383"],
        "urls": [
            ("GSE238219_family.soft.gz", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE238nnn/GSE238219/soft/GSE238219_family.soft.gz", True),
            ("GSE238219_RAW.tar", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE238nnn/GSE238219/suppl/GSE238219_RAW.tar", True),
            ("PMC11865357.html", "https://pmc.ncbi.nlm.nih.gov/articles/PMC11865357/", True),
        ],
    },
    "gse264667": {
        "accessions": ["GSE264667", "PRJNA1100571"],
        "urls": [
            ("GSE264667_family.soft.gz", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE264nnn/GSE264667/soft/GSE264667_family.soft.gz", True),
            ("GSE264667_RAW.tar", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE264nnn/GSE264667/suppl/GSE264667_RAW.tar", False),
        ],
    },
    "gse242934": {
        "accessions": ["GSE242934", "PRJNA1015558"],
        "urls": [
            ("GSE242934_family.soft.gz", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE242nnn/GSE242934/soft/GSE242934_family.soft.gz", True),
            ("GSE242934_Perturbseq_sgrna_lib.txt.gz", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE242nnn/GSE242934/suppl/GSE242934_Perturbseq_sgrna_lib.txt.gz", True),
            ("GSE242934_RAW.tar", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE242nnn/GSE242934/suppl/GSE242934_RAW.tar", False),
        ],
    },
    "cra009621_oscar": {
        "accessions": ["PRJCA014442", "CRA009621"],
        "urls": [
            ("PRJCA014442.html", "https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA014442", True),
            ("CRA009621.html", "https://ngdc.cncb.ac.cn/gsa/browse/CRA009621", True),
            ("OSCAR_article.html", "https://link.springer.com/article/10.1186/s13059-023-03084-8", True),
        ],
    },
    "igvf_hepatocyte": {
        "accessions": ["syn74842722", "Gersbach_WTC11-hepatocyte-differentiation"],
        "urls": [
            ("github_main_commit.json", "https://api.github.com/repos/IGVF/tf_perturb_seq/commits/main", True),
            ("README.md", "https://raw.githubusercontent.com/IGVF/tf_perturb_seq/main/README.md", True),
            ("hepatocyte_README.md", "https://raw.githubusercontent.com/IGVF/tf_perturb_seq/main/datasets/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/README.md", True),
        ],
    },
}


class SourceAuditError(RuntimeError):
    """Raised when a bounded source audit differs or would overwrite evidence."""


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def response_metadata(response: Any) -> dict[str, Any]:
    return {
        "status": getattr(response, "status", None),
        "final_url": response.geturl(),
        "headers": {key.lower(): value for key, value in response.headers.items()},
    }


def acquire(name: str, url: str, target: Path, download: bool) -> dict[str, Any]:
    method = "GET" if download else "HEAD"
    request = Request(url, headers={"User-Agent": "masld-bench-source-audit/1"}, method=method)
    try:
        with urlopen(request, timeout=120) as response:
            result = response_metadata(response)
            result.update({"name": name, "url": url, "method": method})
            if download:
                with target.open("xb") as output:
                    while block := response.read(8 * 1024 * 1024):
                        output.write(block)
                result.update(
                    {
                        "path": target.name,
                        "size_bytes": target.stat().st_size,
                        "sha256": digest(target),
                    }
                )
            return result
    except HTTPError as error:
        return {
            "name": name,
            "url": url,
            "method": method,
            "status": error.code,
            "error": str(error),
            "headers": {key.lower(): value for key, value in error.headers.items()},
        }


def safe_tar_inventory(path: Path) -> list[dict[str, Any]]:
    members: list[dict[str, Any]] = []
    with tarfile.open(path, mode="r:") as archive:
        for member in archive:
            pure = PurePosixPath(member.name)
            if pure.is_absolute() or ".." in pure.parts:
                raise SourceAuditError(f"unsafe archive member: {member.name}")
            members.append(
                {
                    "name": member.name,
                    "size_bytes": member.size,
                    "type": "file" if member.isfile() else "directory" if member.isdir() else "other",
                }
            )
    return members


def audit(source: str, output: Path) -> dict[str, Any]:
    if source not in SOURCES or output.exists():
        raise SourceAuditError("source audit request differs")
    output.mkdir(parents=True, exist_ok=False)
    raw = output / "raw"
    raw.mkdir()
    spec = SOURCES[source]
    records = [
        acquire(name, url, raw / name, download)
        for name, url, download in spec["urls"]
    ]
    inventories: dict[str, Any] = {}
    for record in records:
        path_value = record.get("path")
        if isinstance(path_value, str) and path_value.endswith(".tar"):
            inventories[path_value] = safe_tar_inventory(raw / path_value)
    status = "observed_metadata_only"
    if any(record.get("status") != 200 for record in records):
        status = "observed_with_source_errors"
    result = {
        "schema_version": "masld-bench-open-gene-perturbation-source-audit-v1",
        "source": source,
        "accessions": spec["accessions"],
        "status": status,
        "records": records,
        "archive_inventories": inventories,
        "outcome_model_fit": False,
        "gse313774_accessed": False,
        "inference_warning": "Cells, guides, lanes, GEM groups, wells, and files are not biological replicates.",
    }
    (output / "source_audit.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", choices=tuple(SOURCES), required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.source, args.output)
    print(json.dumps({"source": args.source, "status": result["status"]}, sort_keys=True))


if __name__ == "__main__":
    main()
