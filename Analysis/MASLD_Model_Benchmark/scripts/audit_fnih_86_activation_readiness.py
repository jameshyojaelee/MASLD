#!/usr/bin/env python3
"""Freeze a public-metadata-only readiness audit for the conditional FNIH cohort."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import time
import tomllib
import urllib.error
import urllib.parse
import urllib.request


MAX_SOURCE_BYTES = 32 * 1024 * 1024
VALID_GATE_STATES = {"pass", "fail", "not_evaluable"}
REQUIRED_TOPOLOGIES = {
    "tenx_multiome_rna_atac",
    "paired_tag_rna_h3k27ac",
    "paired_tag_rna_h3k27me3",
    "droplet_hic",
    "visium_hd",
}
ALLOWED_STATIC_ASSET_SUFFIXES = {".js", ".css"}
BIOLOGICAL_FILE_SUFFIXES = {
    ".bed",
    ".bedgraph",
    ".bigwig",
    ".bw",
    ".cool",
    ".csv",
    ".feather",
    ".gz",
    ".h5",
    ".h5ad",
    ".hic",
    ".loom",
    ".mtx",
    ".parquet",
    ".rds",
    ".tsv",
    ".vcf",
    ".zarr",
}
PORTAL_ASSET_PATTERN = re.compile(
    r"(?:src|href)=[\"']([^\"']+\.(?:js|css)(?:\?[^\"']*)?)[\"']",
    re.IGNORECASE,
)
QUOTED_PATH_PATTERN = re.compile(
    r"[\"']([^\"'\s]{1,500}\.(?:bed|bedgraph|bigwig|bw|cool|csv|feather|gz|h5|h5ad|hic|loom|mtx|parquet|rds|tsv|vcf|zarr)(?:\?[^\"']*)?)[\"']",
    re.IGNORECASE,
)


class FNIHReadinessError(RuntimeError):
    """Raised when the readiness requirements or public evidence is invalid."""


def canonical_json(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def validate_contract(contract: dict[str, object]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    for key in (
        "automatic_activation",
        "sealed_outcomes_read",
        "biological_matrices_downloaded",
        "controlled_data_accessed",
    ):
        if contract.get(key) is not False:
            raise FNIHReadinessError(f"{key} must remain false")
    if contract.get("expected_final_donor_n") != 86:
        raise FNIHReadinessError("final donor count must remain 86")
    if contract.get("source_collected_donor_n") != 87:
        raise FNIHReadinessError("source-collected donor count must remain 87")
    if contract.get("source_excluded_qc_donor_n") != 1:
        raise FNIHReadinessError("one source donor must remain recorded as excluded")
    publication = contract.get("publication")
    if not isinstance(publication, dict) or publication.get("pmcid") != "PMC12083587":
        raise FNIHReadinessError("publication authority is missing or changed")
    sources = contract.get("public_source")
    if not isinstance(sources, list) or len(sources) < 7:
        raise FNIHReadinessError("public evidence source census is incomplete")
    source_ids = [row.get("source_id") for row in sources]
    if len(source_ids) != len(set(source_ids)):
        raise FNIHReadinessError("public source IDs must be unique")
    for row in sources:
        url = str(row.get("url", ""))
        parsed = urllib.parse.urlsplit(url)
        if parsed.scheme != "https" or not parsed.netloc:
            raise FNIHReadinessError(f"public source must be HTTPS: {url}")
        if not row.get("required_tokens"):
            raise FNIHReadinessError(f"required evidence tokens missing: {row.get('source_id')}")
    topologies = contract.get("assay_topology")
    if not isinstance(topologies, list):
        raise FNIHReadinessError("assay topology census is missing")
    topology_ids = {str(row.get("assay_id")) for row in topologies}
    if topology_ids != REQUIRED_TOPOLOGIES:
        raise FNIHReadinessError(f"assay topology mismatch: {sorted(topology_ids)}")
    by_id = {str(row["assay_id"]): row for row in topologies}
    if by_id["tenx_multiome_rna_atac"].get("within_assay_pairing") != "same_nucleus":
        raise FNIHReadinessError("10x Multiome must remain same-nucleus RNA-ATAC")
    if by_id["tenx_multiome_rna_atac"].get("donor_n") != 86:
        raise FNIHReadinessError("10x Multiome donor count must remain 86")
    if by_id["droplet_hic"].get("donor_n") != 85:
        raise FNIHReadinessError("droplet Hi-C donor count must remain 85")
    if by_id["visium_hd"].get("donor_n") != 7:
        raise FNIHReadinessError("Visium HD donor count must remain 7")
    if any(row.get("public_donor_barcode_join") is not False for row in topologies):
        raise FNIHReadinessError("no public donor-barcode join may be asserted")
    criteria = contract.get("activation_criterion")
    if not isinstance(criteria, list) or not criteria:
        raise FNIHReadinessError("activation criteria are missing")
    criterion_ids = [row.get("criterion_id") for row in criteria]
    if len(criterion_ids) != len(set(criterion_ids)):
        raise FNIHReadinessError("activation criterion IDs must be unique")
    for row in criteria:
        if row.get("required") is not True:
            raise FNIHReadinessError("every v1 activation criterion must be required")
        if row.get("current_state") not in VALID_GATE_STATES:
            raise FNIHReadinessError(f"invalid gate state: {row.get('criterion_id')}")
    decision = contract.get("decision")
    if not isinstance(decision, dict):
        raise FNIHReadinessError("decision block is missing")
    if decision.get("current_activation_state") != "blocked":
        raise FNIHReadinessError("FNIH source must remain blocked in this contract")
    if decision.get("eligible_as_external_seal_now") is not False:
        raise FNIHReadinessError("FNIH cannot be externally seal-eligible now")
    if "gse289173" in json.dumps(contract).lower():
        raise FNIHReadinessError("sealed accession must not enter the FNIH audit")
    return sources, criteria


def evaluate_gate(criteria: list[dict[str, object]]) -> dict[str, object]:
    required = [row for row in criteria if row.get("required") is True]
    nonpassing = [
        {
            "criterion_id": row["criterion_id"],
            "state": row["current_state"],
            "evidence": row["evidence"],
        }
        for row in required
        if row.get("current_state") != "pass"
    ]
    return {
        "schema_version": "masld-bench-fnih-activation-gate-v1",
        "automatic_activation": False,
        "required_criterion_count": len(required),
        "passing_criterion_count": len(required) - len(nonpassing),
        "nonpassing_criterion_count": len(nonpassing),
        "nonpassing": nonpassing,
        "activation_gate_pass": not nonpassing,
        "current_activation_state": "eligible_for_manual_manifest_revision" if not nonpassing else "blocked",
        "eligible_as_external_seal_now": False,
        "sealed_outcomes_read": False,
        "biological_matrices_downloaded": False,
        "controlled_data_accessed": False,
    }


def _safe_static_asset_url(portal_url: str, candidate: str) -> str:
    absolute = urllib.parse.urljoin(portal_url, candidate)
    portal = urllib.parse.urlsplit(portal_url)
    parsed = urllib.parse.urlsplit(absolute)
    suffix = Path(parsed.path).suffix.lower()
    if parsed.scheme != "https" or parsed.netloc != portal.netloc:
        raise FNIHReadinessError(f"portal asset left the source origin: {absolute}")
    if not parsed.path.startswith("/MASLD/assets/") or suffix not in ALLOWED_STATIC_ASSET_SUFFIXES:
        raise FNIHReadinessError(f"non-code portal asset request rejected: {absolute}")
    return absolute


def discover_portal_assets(portal_url: str, html: str) -> list[str]:
    return sorted({_safe_static_asset_url(portal_url, path) for path in PORTAL_ASSET_PATTERN.findall(html)})


def inventory_candidate_data_paths(asset_id: str, payload: bytes) -> list[dict[str, object]]:
    text = payload.decode("utf-8", errors="replace")
    rows: list[dict[str, object]] = []
    for candidate in sorted(set(QUOTED_PATH_PATTERN.findall(text))):
        parsed = urllib.parse.urlsplit(candidate)
        suffix = Path(parsed.path).suffix.lower()
        if suffix not in BIOLOGICAL_FILE_SUFFIXES:
            continue
        rows.append(
            {
                "source_asset_id": asset_id,
                "literal": candidate,
                "suffix": suffix,
                "requested": False,
                "interpretation": "public_code_string_only_not_a_verified_file_or_rights_manifest",
            }
        )
    return rows


def _read_response(response: object, source_id: str) -> tuple[bytes, dict[str, object]]:
    payload = response.read(MAX_SOURCE_BYTES + 1)
    if len(payload) > MAX_SOURCE_BYTES:
        raise FNIHReadinessError(f"public metadata source exceeds cap: {source_id}")
    headers = response.headers
    return payload, {
        "source_id": source_id,
        "url": response.geturl(),
        "http_status": response.status,
        "retrieved_at_utc": utc_now(),
        "http_date": headers.get("Date", ""),
        "content_type": headers.get("Content-Type", ""),
        "last_modified": headers.get("Last-Modified", ""),
        "etag": headers.get("ETag", ""),
        "bytes": len(payload),
        "sha256": sha256_bytes(payload),
    }


def fetch_public_source(source_id: str, url: str) -> tuple[bytes, dict[str, object]]:
    parsed = urllib.parse.urlsplit(url)
    if parsed.scheme != "https":
        raise FNIHReadinessError(f"non-HTTPS request rejected: {url}")
    request = urllib.request.Request(
        url,
        headers={
            "Accept": "application/json,application/xml,text/html,text/plain,text/css,application/javascript,*/*",
            "User-Agent": "masld-bench-fnih-public-metadata-audit/1.0",
        },
    )
    last_error: Exception | None = None
    for attempt in range(1, 4):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                payload, receipt = _read_response(response, source_id)
                receipt["attempt"] = attempt
                return payload, receipt
        except (urllib.error.URLError, TimeoutError, ConnectionError) as error:
            last_error = error
            if attempt < 3:
                time.sleep(2 ** (attempt - 1))
    raise FNIHReadinessError(f"failed to fetch {source_id}: {last_error}") from last_error


def require_tokens(source_id: str, payload: bytes, tokens: list[object]) -> None:
    text = payload.decode("utf-8", errors="replace").lower()
    missing = [str(token) for token in tokens if str(token).lower() not in text]
    if missing:
        raise FNIHReadinessError(f"{source_id} missing public evidence tokens: {missing}")


def safe_filename(source_id: str, url: str) -> str:
    suffix = Path(urllib.parse.urlsplit(url).path).suffix.lower()
    if suffix not in {".css", ".html", ".js", ".json", ".xml"}:
        suffix = ".txt"
    return re.sub(r"[^a-z0-9_.-]+", "_", source_id.lower()) + suffix


def fetch_github_inventory(
    repository_payload: bytes,
) -> list[tuple[str, str]]:
    record = json.loads(repository_payload)
    full_name = record.get("full_name")
    default_branch = record.get("default_branch")
    if full_name != "Gaulton-Lab/FNIH.Liver" or not default_branch:
        raise FNIHReadinessError("GitHub repository identity or default branch is unresolved")
    quoted_branch = urllib.parse.quote(str(default_branch), safe="")
    commit_url = f"https://api.github.com/repos/{full_name}/commits/{quoted_branch}"
    return [("github_default_branch_commit", commit_url)]


def fetch_github_tree(commit_payload: bytes) -> tuple[str, str]:
    record = json.loads(commit_payload)
    commit_sha = str(record.get("sha", ""))
    if not re.fullmatch(r"[0-9a-f]{40}", commit_sha):
        raise FNIHReadinessError("GitHub default-branch commit SHA is invalid")
    return commit_sha, f"https://api.github.com/repos/Gaulton-Lab/FNIH.Liver/git/trees/{commit_sha}?recursive=1"


def summarize_public_evidence(
    receipts: list[dict[str, object]],
    payloads: dict[str, bytes],
    candidate_paths: list[dict[str, object]],
) -> dict[str, object]:
    zenodo = json.loads(payloads["zenodo_record"])
    zenodo_metadata = zenodo.get("metadata", {})
    zenodo_files = zenodo.get("files") or []
    github = json.loads(payloads["github_repository"]) if "github_repository" in payloads else {}
    github_commit = (
        json.loads(payloads["github_default_branch_commit"])
        if "github_default_branch_commit" in payloads
        else {}
    )
    github_tree = (
        json.loads(payloads["github_repository_tree"])
        if "github_repository_tree" in payloads
        else {}
    )
    gds_results: dict[str, object] = {}
    for source_id in ("gds_by_pmid", "gds_by_exact_title"):
        result = json.loads(payloads[source_id]).get("esearchresult", {})
        gds_results[source_id] = {
            "count": int(result.get("count", 0)),
            "ids": list(result.get("idlist", [])),
            "interpretation": "candidate_repository_matches_only_no_accession_is_auto_admitted",
        }
    license_record = zenodo_metadata.get("license") or {}
    suffix_counts: dict[str, int] = {}
    modality_counts: dict[str, int] = {}
    identity_keywords = ("donor", "participant", "sample", "barcode")
    for row in candidate_paths:
        suffix = str(row["suffix"])
        suffix_counts[suffix] = suffix_counts.get(suffix, 0) + 1
        literal = str(row["literal"])
        match = re.search(r"/FNIH\.Liver/(ATAC|H3K27ac|H3K27me3)/", literal)
        if match:
            modality = match.group(1)
            modality_counts[modality] = modality_counts.get(modality, 0) + 1
    return {
        "schema_version": "masld-bench-fnih-public-evidence-summary-v1",
        "source_receipt_count": len(receipts),
        "all_requests_public_metadata_or_static_code": True,
        "portal_static_asset_count": sum(str(row["source_id"]).startswith("portal_asset_") for row in receipts),
        "portal_candidate_data_path_string_count": len(candidate_paths),
        "portal_candidate_data_paths_requested": False,
        "portal_candidate_suffix_counts": suffix_counts,
        "portal_candidate_modality_counts": modality_counts,
        "portal_aggregate_celltype_or_disease_path_count": sum(
            "/celltypes/" in str(row["literal"]) or "/disease/" in str(row["literal"])
            for row in candidate_paths
        ),
        "portal_identity_keyword_path_count": sum(
            any(token in str(row["literal"]).lower() for token in identity_keywords)
            for row in candidate_paths
        ),
        "zenodo": {
            "record_id": zenodo.get("id"),
            "access_right": zenodo_metadata.get("access_right"),
            "license_id": license_record.get("id"),
            "publication_date": zenodo_metadata.get("publication_date"),
            "file_count": len(zenodo_files),
            "rights_interpretation": "record_metadata_and_license_do_not_establish_rights_to_restricted_or_absent_biological_files",
        },
        "github": {
            "metadata_available": bool(github),
            "full_name": github.get("full_name"),
            "default_branch": github.get("default_branch"),
            "default_branch_commit_sha": github_commit.get("sha"),
            "repository_license": (github.get("license") or {}).get("spdx_id"),
            "tree_entry_count": len(github_tree.get("tree") or []) if github_tree else None,
            "tree_truncated": github_tree.get("truncated"),
            "rights_interpretation": "code_repository_terms_are_separate_from_participant_data_rights",
        },
        "gds_search": gds_results,
        "sealed_outcomes_read": False,
        "biological_matrices_downloaded": False,
        "controlled_data_accessed": False,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    contract_path = args.contract.resolve(strict=True)
    output = args.output.resolve()
    if output.exists():
        raise FNIHReadinessError(f"refusing to overwrite output: {output}")
    contract = tomllib.loads(contract_path.read_text(encoding="utf-8"))
    sources, criteria = validate_contract(contract)
    metadata_dir = output / "public_metadata"
    metadata_dir.mkdir(parents=True)

    receipts: list[dict[str, object]] = []
    payloads: dict[str, bytes] = {}
    for spec in sources:
        source_id = str(spec["source_id"])
        url = str(spec["url"])
        try:
            payload, receipt = fetch_public_source(source_id, url)
        except FNIHReadinessError as error:
            if spec.get("allow_unavailable") is not True:
                raise
            receipts.append(
                {
                    "source_id": source_id,
                    "url": url,
                    "content_class": spec["content_class"],
                    "available": False,
                    "retrieved_at_utc": utc_now(),
                    "error_class": type(error.__cause__).__name__ if error.__cause__ else type(error).__name__,
                    "error": str(error),
                    "unavailable_disposition": spec["unavailable_disposition"],
                }
            )
            continue
        require_tokens(source_id, payload, list(spec["required_tokens"]))
        receipt["available"] = True
        receipt["content_class"] = spec["content_class"]
        receipt["required_evidence_pass"] = True
        target = metadata_dir / safe_filename(source_id, url)
        target.write_bytes(payload)
        receipt["local_path"] = str(target.relative_to(output))
        receipts.append(receipt)
        payloads[source_id] = payload

    portal_url = next(str(row["url"]) for row in sources if row["source_id"] == "portal_html")
    portal_assets = discover_portal_assets(
        portal_url, payloads["portal_html"].decode("utf-8", errors="replace")
    )
    candidate_paths: list[dict[str, object]] = []
    portal_asset_inventory: list[dict[str, object]] = []
    for index, asset_url in enumerate(portal_assets, start=1):
        source_id = f"portal_asset_{index:02d}"
        payload, receipt = fetch_public_source(source_id, asset_url)
        receipt["content_class"] = "public_static_code_asset"
        receipt["available"] = True
        target = metadata_dir / safe_filename(source_id, asset_url)
        target.write_bytes(payload)
        receipt["local_path"] = str(target.relative_to(output))
        receipts.append(receipt)
        payloads[source_id] = payload
        asset_candidates = inventory_candidate_data_paths(source_id, payload)
        candidate_paths.extend(asset_candidates)
        portal_asset_inventory.append(
            {
                "source_id": source_id,
                "url": asset_url,
                "sha256": receipt["sha256"],
                "bytes": receipt["bytes"],
                "candidate_data_path_string_count": len(asset_candidates),
                "biological_file_requested": False,
            }
        )

    if "github_repository" in payloads:
        try:
            for source_id, url in fetch_github_inventory(payloads["github_repository"]):
                payload, receipt = fetch_public_source(source_id, url)
                receipt["content_class"] = "public_code_metadata_json"
                receipt["available"] = True
                target = metadata_dir / safe_filename(source_id, url)
                target.write_bytes(payload)
                receipt["local_path"] = str(target.relative_to(output))
                receipts.append(receipt)
                payloads[source_id] = payload
            _, tree_url = fetch_github_tree(payloads["github_default_branch_commit"])
            tree_payload, tree_receipt = fetch_public_source("github_repository_tree", tree_url)
            tree_receipt["content_class"] = "public_code_tree_metadata_json"
            tree_receipt["available"] = True
            tree_target = metadata_dir / "github_repository_tree.json"
            tree_target.write_bytes(tree_payload)
            tree_receipt["local_path"] = str(tree_target.relative_to(output))
            receipts.append(tree_receipt)
            payloads["github_repository_tree"] = tree_payload
        except FNIHReadinessError as error:
            receipts.append(
                {
                    "source_id": "github_dynamic_inventory",
                    "url": "https://api.github.com/repos/Gaulton-Lab/FNIH.Liver",
                    "content_class": "public_code_metadata_json",
                    "available": False,
                    "retrieved_at_utc": utc_now(),
                    "error_class": type(error.__cause__).__name__ if error.__cause__ else type(error).__name__,
                    "error": str(error),
                    "unavailable_disposition": "Exact code commit or tree remains unresolved; participant-data activation stays blocked.",
                }
            )

    gate = evaluate_gate(criteria)
    public_summary = summarize_public_evidence(receipts, payloads, candidate_paths)
    topology_rights_reference = {
        "schema_version": "masld-bench-fnih-topology-rights-reference-v1",
        "dataset_id": contract["dataset_id"],
        "expected_final_donor_n": contract["expected_final_donor_n"],
        "source_collected_donor_n": contract["source_collected_donor_n"],
        "source_excluded_qc_donor_n": contract["source_excluded_qc_donor_n"],
        "assay_topology": contract["assay_topology"],
        "native_reference_high_level": contract["native_reference_high_level"],
        "native_reference_exact_bundle": contract["native_reference_exact_bundle"],
        "native_annotation_release": contract["native_annotation_release"],
        "rights_interpretation": contract["rights_interpretation"],
        "same_nucleus_scope": "within Multiome and within each Paired-Tag mark only; no cross-assay nucleus pairing",
        "public_donor_barcode_join_available": False,
        "exact_reference_ready": False,
    }
    (output / "source_receipts.json").write_bytes(canonical_json({"sources": receipts}))
    (output / "portal_asset_inventory.json").write_bytes(
        canonical_json(
            {
                "assets": portal_asset_inventory,
                "candidate_data_path_strings": candidate_paths,
                "candidate_paths_requested": False,
            }
        )
    )
    (output / "public_evidence_summary.json").write_bytes(canonical_json(public_summary))
    (output / "topology_reference_rights.json").write_bytes(
        canonical_json(topology_rights_reference)
    )
    (output / "fnih_activation_gate.json").write_bytes(canonical_json(gate))
    (output / "DECISION.md").write_text(
        "# FNIH 86-liver activation-readiness decision\n\n"
        "**Decision: blocked.** The public paper establishes a high-value same-nucleus "
        "RNA-ATAC topology for 86 retained donors, but the public portal is not by itself "
        "an admitted donor-resolved dataset.\n\n"
        f"The deterministic gate has {gate['passing_criterion_count']} passing and "
        f"{gate['nonpassing_criterion_count']} nonpassing required criteria. Activation "
        "requires every criterion to pass in one immutable audit, followed by a reviewed "
        "DatasetManifest and TaskSpec revision. No biological matrix, controlled data, or "
        "sealed outcome was accessed.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
