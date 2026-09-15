#!/usr/bin/env python3
"""Exact inclusion requirements for the frozen GSE105127 legacy reference."""

from __future__ import annotations

from typing import Any, Mapping


REFERENCE_ARTIFACT_CLASS = "gse105127_reference_bundle"
REFERENCE_RECEIPT_SCHEMA = "masld-bench-gse105127-reference-bundle-v2"
REFERENCE_RECEIPT_STATUS = "passed_exact_legacy_stream"
REFERENCE_SOURCE_ASSEMBLY = "1000_Genomes_GRCh37_human_g1k_v37"


def admits_exact_reference(
    manifest: Mapping[str, Any], receipt: Mapping[str, Any]
) -> bool:
    metadata = manifest.get("metadata", {})
    return bool(
        metadata.get("artifact_class") == REFERENCE_ARTIFACT_CLASS
        and metadata.get("status") == REFERENCE_RECEIPT_STATUS
        and receipt.get("schema_version") == REFERENCE_RECEIPT_SCHEMA
        and receipt.get("status") == REFERENCE_RECEIPT_STATUS
        and receipt.get("source_assembly") == REFERENCE_SOURCE_ASSEMBLY
        and receipt.get("labels_accessed") is False
        and receipt.get("fit_or_score_performed") is False
    )
