#!/usr/bin/env python3
"""Independent validator for an emitted MASLD Gene Catalog v2 contract."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

import gene_catalog_v2 as catalog


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate(candidate_dir: Path, review_dir: Path) -> None:
    if review_dir.exists():
        raise FileExistsError(f"refusing to overwrite {review_dir}")
    review_dir.mkdir(parents=True)

    manifest = json.loads((candidate_dir / "manifest.json").read_text())
    assert manifest["catalog_schema_version"] == catalog.CATALOG_SCHEMA_VERSION
    assert manifest["contains_gene_outcomes"] is False
    assert manifest["genetics_dependent_entries"] == "blocked_pending_corrected_coloc"
    for artifact in manifest["artifacts"]:
        path = candidate_dir / artifact["path"]
        assert path.stat().st_size == artifact["bytes"]
        assert _sha256(path) == artifact["sha256"]

    schema = json.loads(
        (candidate_dir / "masld_gene_catalog_v2_schema.json").read_text()
    )
    dictionary = pd.read_csv(
        candidate_dir / "masld_gene_catalog_v2_data_dictionary.tsv", sep="\t"
    )
    rules = pd.read_csv(
        candidate_dir / "masld_gene_catalog_v2_experiment_rules.tsv", sep="\t"
    )
    assert set(catalog.V2_FIELDS).issubset(set(dictionary["field"]))
    assert set(dictionary.loc[dictionary["required_v2"], "field"]) == set(
        schema["required"]
    )
    assert set(rules["experiment_rule_id"]) == set(catalog.EXPERIMENT_RULES_V2)
    assert not any("screen" in column.lower() for column in dictionary["field"])
    assert "top_snp" not in "\n".join(dictionary.astype(str).stack()).lower()
    assert "no gene outcomes" in (candidate_dir / "README.md").read_text().lower()

    verdict = {
        "status": "pass",
        "candidate_dir": str(candidate_dir.resolve()),
        "catalog_schema_version": catalog.CATALOG_SCHEMA_VERSION,
        "field_count": int(len(dictionary)),
        "rule_count": int(len(rules)),
        "contains_gene_outcomes": False,
        "genetics_dependent_entries": "blocked_pending_corrected_coloc",
        "candidate_manifest_sha256": _sha256(candidate_dir / "manifest.json"),
    }
    with (review_dir / "VALIDATED.json").open("w", encoding="utf-8") as handle:
        json.dump(verdict, handle, indent=2, sort_keys=True)
        handle.write("\n")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-dir", type=Path, required=True)
    parser.add_argument("--review-dir", type=Path, required=True)
    args = parser.parse_args()
    validate(args.candidate_dir.resolve(), args.review_dir.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
