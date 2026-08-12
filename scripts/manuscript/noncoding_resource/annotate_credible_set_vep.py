#!/usr/bin/env python3
"""Annotate the complete unique credible-set variant universe with GRCh37 VEP."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import time
import urllib.error
import urllib.request

SERVER = "https://grch37.rest.ensembl.org"
ENDPOINT = "/vep/homo_sapiens/region"
MAX_BATCH = 200
MAX_ATTEMPTS = 8

SEVERITY = {
    consequence: rank
    for rank, consequence in enumerate(
        (
            "transcript_ablation",
            "splice_acceptor_variant",
            "splice_donor_variant",
            "stop_gained",
            "frameshift_variant",
            "stop_lost",
            "start_lost",
            "transcript_amplification",
            "inframe_insertion",
            "inframe_deletion",
            "missense_variant",
            "protein_altering_variant",
            "splice_region_variant",
            "incomplete_terminal_codon_variant",
            "start_retained_variant",
            "stop_retained_variant",
            "synonymous_variant",
            "coding_sequence_variant",
            "mature_miRNA_variant",
            "5_prime_UTR_variant",
            "3_prime_UTR_variant",
            "non_coding_transcript_exon_variant",
            "intron_variant",
            "NMD_transcript_variant",
            "non_coding_transcript_variant",
            "upstream_gene_variant",
            "downstream_gene_variant",
            "TFBS_ablation",
            "TFBS_amplification",
            "TF_binding_site_variant",
            "regulatory_region_ablation",
            "regulatory_region_amplification",
            "feature_elongation",
            "regulatory_region_variant",
            "feature_truncation",
            "intergenic_variant",
        )
    )
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--members", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--batch-size", type=int, default=MAX_BATCH)
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_variants(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"variant_id", "chromosome", "position", "allele1", "allele2"}
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise RuntimeError(f"missing member fields: {sorted(missing)}")
        by_id: dict[str, dict[str, str]] = {}
        for row in reader:
            compact = {field: row[field] for field in required}
            previous = by_id.setdefault(row["variant_id"], compact)
            if previous != compact:
                raise RuntimeError(f"variant identity drift: {row['variant_id']}")
    return [by_id[key] for key in sorted(by_id)]


def vcf_line(row: dict[str, str]) -> str:
    # Upstream fine-mapping convention is allele1=effect/ALT, allele2=other/REF.
    return " ".join(
        (
            row["chromosome"],
            row["position"],
            row["variant_id"],
            row["allele2"],
            row["allele1"],
            ".",
            ".",
            ".",
        )
    )


def post_batch(variants: list[str]) -> tuple[list[dict[str, object]], str]:
    body = json.dumps({"variants": variants}, separators=(",", ":")).encode()
    request = urllib.request.Request(
        SERVER + ENDPOINT,
        data=body,
        headers={"Content-Type": "application/json", "Accept": "application/json"},
        method="POST",
    )
    for attempt in range(MAX_ATTEMPTS):
        try:
            with urllib.request.urlopen(request, timeout=180) as response:
                raw = response.read()
                return json.loads(raw), response.headers.get("X-RateLimit-Remaining", "")
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as error:
            if attempt + 1 == MAX_ATTEMPTS:
                raise RuntimeError(f"VEP request failed after {MAX_ATTEMPTS} attempts") from error
            retry_after = getattr(error, "headers", {}).get("Retry-After") if hasattr(error, "headers") else None
            delay = float(retry_after) if retry_after else min(60.0, 2.0**attempt)
            time.sleep(delay)
    raise AssertionError("unreachable")


def classify(consequences: set[str]) -> tuple[str, str]:
    ordered = sorted(consequences, key=lambda term: (SEVERITY.get(term, 10_000), term))
    if not ordered:
        return "", ""
    if any(
        term
        in {
            "transcript_ablation",
            "stop_gained",
            "frameshift_variant",
            "stop_lost",
            "start_lost",
            "transcript_amplification",
            "inframe_insertion",
            "inframe_deletion",
            "missense_variant",
            "protein_altering_variant",
        }
        for term in consequences
    ):
        source_class = "coding_protein_altering"
    elif "synonymous_variant" in consequences:
        source_class = "coding_synonymous"
    elif any("splice" in term for term in consequences):
        source_class = "splice_region"
    else:
        source_class = "noncoding"
    return ",".join(ordered), source_class


def write_tsv(path: Path, rows: list[dict[str, object]], fields: tuple[str, ...]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    if args.output_dir.exists():
        raise RuntimeError(f"output exists: {args.output_dir}")
    if not 1 <= args.batch_size <= 300:
        raise RuntimeError("batch size must be in [1,300]")
    variants = read_variants(args.members)
    temporary = args.output_dir.with_name(
        f".{args.output_dir.name}.tmp.{os.environ.get('SLURM_JOB_ID', os.getpid())}"
    )
    temporary.mkdir(parents=True, exist_ok=False)
    try:
        records: list[dict[str, object]] = []
        batch_manifest: list[dict[str, object]] = []
        for batch_index, start in enumerate(range(0, len(variants), args.batch_size), 1):
            batch = variants[start : start + args.batch_size]
            submitted = [vcf_line(row) for row in batch]
            response, remaining = post_batch(submitted)
            by_id = {str(item.get("id", "")): item for item in response}
            if set(by_id) != {row["variant_id"] for row in batch}:
                missing = sorted({row["variant_id"] for row in batch} - set(by_id))
                extra = sorted(set(by_id) - {row["variant_id"] for row in batch})
                raise RuntimeError(f"VEP response identity mismatch; missing={missing[:3]} extra={extra[:3]}")
            for row in batch:
                item = by_id[row["variant_id"]]
                transcript_terms = {
                    term
                    for consequence in item.get("transcript_consequences", [])
                    for term in consequence.get("consequence_terms", [])
                }
                regulatory_terms = {
                    term
                    for consequence in item.get("regulatory_feature_consequences", [])
                    for term in consequence.get("consequence_terms", [])
                }
                intergenic_terms = set(item.get("most_severe_consequence", "").split("&"))
                consequences = transcript_terms | regulatory_terms | intergenic_terms
                consequences.discard("")
                consequence, source_class = classify(consequences)
                records.append(
                    {
                        "variant_id": row["variant_id"],
                        "chromosome": row["chromosome"],
                        "position": row["position"],
                        "effect_allele": row["allele1"],
                        "other_allele": row["allele2"],
                        "vep_allele_string": item.get("allele_string", ""),
                        "most_severe_consequence": item.get("most_severe_consequence", ""),
                        "Consequence": consequence,
                        "class": source_class,
                        "annotation_source": "Ensembl_GRCh37_REST_VEP",
                    }
                )
            batch_manifest.append(
                {
                    "batch": batch_index,
                    "n_submitted": len(batch),
                    "n_returned": len(response),
                    "rate_limit_remaining": remaining,
                }
            )
            if batch_index % 10 == 0:
                print(f"VEP_PROGRESS\t{min(start + args.batch_size, len(variants))}/{len(variants)}", flush=True)

        records.sort(key=lambda row: str(row["variant_id"]))
        write_tsv(
            temporary / "credible_set_variant_consequences_complete.tsv",
            records,
            (
                "variant_id",
                "chromosome",
                "position",
                "effect_allele",
                "other_allele",
                "vep_allele_string",
                "most_severe_consequence",
                "Consequence",
                "class",
                "annotation_source",
            ),
        )
        write_tsv(
            temporary / "request_manifest.tsv",
            batch_manifest,
            ("batch", "n_submitted", "n_returned", "rate_limit_remaining"),
        )
        write_tsv(
            temporary / "source_manifest.tsv",
            [
                {
                    "source_path": str(args.members.resolve()),
                    "size_bytes": args.members.stat().st_size,
                    "sha256": sha256_file(args.members),
                    "server": SERVER,
                    "endpoint": ENDPOINT,
                    "n_unique_variants": len(variants),
                    "allele_convention": "allele2=REF;allele1=ALT",
                }
            ],
            (
                "source_path",
                "size_bytes",
                "sha256",
                "server",
                "endpoint",
                "n_unique_variants",
                "allele_convention",
            ),
        )
        with (temporary / "execution_manifest.json").open("x", encoding="utf-8") as handle:
            json.dump(
                {
                    "producer": str(Path(__file__).resolve()),
                    "producer_sha256": sha256_file(Path(__file__).resolve()),
                    "n_unique_variants": len(variants),
                    "n_annotated_variants": len(records),
                    "batch_size": args.batch_size,
                    "server": SERVER,
                    "endpoint": ENDPOINT,
                },
                handle,
                indent=2,
                sort_keys=True,
            )
            handle.write("\n")
        outputs = []
        for path in sorted(temporary.iterdir()):
            outputs.append(
                {
                    "relative_path": path.name,
                    "size_bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                }
            )
        write_tsv(
            temporary / "output_manifest.tsv",
            outputs,
            ("relative_path", "size_bytes", "sha256"),
        )
        args.output_dir.parent.mkdir(parents=True, exist_ok=True)
        temporary.rename(args.output_dir)
    except Exception:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    print(json.dumps({"output": str(args.output_dir), "n_variants": len(records)}, sort_keys=True))


if __name__ == "__main__":
    main()
