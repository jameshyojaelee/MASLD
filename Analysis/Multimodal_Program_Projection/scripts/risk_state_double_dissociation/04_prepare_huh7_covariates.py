#!/usr/bin/env python3
"""Extract outcome-independent HuH-7 observability covariates for CLCC1 matching."""

from __future__ import annotations

import csv
import gzip
import json
import re
from pathlib import Path

from risk_state_common import (
    CANDIDATE_ROOT,
    PROJECT_ROOT,
    atomic_write_text,
    require_validated_seal,
    sha256_file,
    write_tsv,
)


DEPMAP_ROOT = (
    PROJECT_ROOT
    / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/myojin_hlf/source/depmap_24q4"
)
EXPECTED = {
    "Model.csv": "b7a0c1385e6cef30132b56aff61f1261d11e3f490490b355c430d32ee0dbdcfa",
    "CRISPRGeneEffect.csv": "3d8f3ec6dbf2db7ff834b79b508622ec0b226f3518003fe96ecf5a4fcf167e3b",
    "OmicsExpressionProteinCodingGenesTPMLogp1.csv": "2a71dc94110efcc0221eae821bb93a9f03b54bea16f005818911a09d33383d56",
}
MODEL_ID = "ACH-000480"


def parse_gene_header(value: str) -> tuple[str, str]:
    match = re.match(r"^(.*) \(([0-9]+)\)$", value)
    if match:
        return match.group(1), match.group(2)
    return value, ""


def extract_model_row(path: Path) -> dict[tuple[str, str], str]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        for row in reader:
            if row and row[0] == MODEL_ID:
                if len(row) != len(header):
                    raise RuntimeError(f"Column mismatch in {path.name}: {len(row)} != {len(header)}")
                return {parse_gene_header(column): value for column, value in zip(header[1:], row[1:])}
    raise RuntimeError(f"HuH-7 model {MODEL_ID} absent from {path}")


def gene_annotations() -> tuple[dict[str, str], dict[str, int]]:
    biotypes: dict[str, str] = {}
    with gzip.open(PROJECT_ROOT / "data/gencode_v49_gene_metadata.tsv.gz", "rt", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            biotypes.setdefault(row["gene_name"], row["gene_biotype"])

    lengths: dict[str, int] = {}
    gtf = Path("/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz")
    with gzip.open(gtf, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            match = re.search(r'gene_name "([^"]+)"', fields[8])
            if match:
                length = int(fields[4]) - int(fields[3]) + 1
                lengths[match.group(1)] = max(lengths.get(match.group(1), 0), length)
    return biotypes, lengths


def numeric(value: str) -> float | str:
    try:
        return float(value)
    except (TypeError, ValueError):
        return ""


def main() -> None:
    seal = require_validated_seal()
    observed: dict[str, str] = {}
    for filename, expected in EXPECTED.items():
        path = DEPMAP_ROOT / filename
        digest = sha256_file(path)
        if digest != expected:
            raise RuntimeError(f"DepMap source drift for {filename}: {digest} != {expected}")
        observed[filename] = digest

    expression = extract_model_row(DEPMAP_ROOT / "OmicsExpressionProteinCodingGenesTPMLogp1.csv")
    chronos = extract_model_row(DEPMAP_ROOT / "CRISPRGeneEffect.csv")
    biotypes, lengths = gene_annotations()
    keys = sorted(set(expression) | set(chronos))
    rows = []
    for symbol, entrez in keys:
        log2_tpm1 = numeric(expression.get((symbol, entrez), ""))
        tpm = 2**log2_tpm1 - 1 if isinstance(log2_tpm1, float) else ""
        rows.append(
            {
                "gene_symbol": symbol,
                "entrez_id": entrez,
                "huh7_expression_log2_tpm_plus_1": log2_tpm1,
                "huh7_tpm": tpm,
                "huh7_chronos": numeric(chronos.get((symbol, entrez), "")),
                "gencode_v49_gene_length": lengths.get(symbol, ""),
                "gencode_v49_biotype": biotypes.get(symbol, ""),
                "depmap_model_id": MODEL_ID,
                "depmap_release": "DepMap Public 24Q4",
            }
        )
    output = CANDIDATE_ROOT / "preprocessed/huh7_observability_covariates.tsv"
    write_tsv(
        output,
        rows,
        [
            "gene_symbol", "entrez_id", "huh7_expression_log2_tpm_plus_1", "huh7_tpm",
            "huh7_chronos", "gencode_v49_gene_length", "gencode_v49_biotype",
            "depmap_model_id", "depmap_release",
        ],
    )
    amendment = {
        "amendment_type": "outcome_independent_covariate_source_identity",
        "reason": "Plan 42 froze HuH-7 expression and generic fitness as matching covariates; their checksum-pinned predecessor objects are referenced rather than copied",
        "external_outcomes_used_to_choose_source": False,
        "model_id": MODEL_ID,
        "source_sha256": observed,
        "output_sha256": sha256_file(output),
        "specification_sha256": seal["specification_sha256"],
    }
    atomic_write_text(
        CANDIDATE_ROOT / "COVARIATE_SOURCE_AMENDMENT.json",
        json.dumps(amendment, indent=2, sort_keys=True) + "\n",
    )
    print(json.dumps({**amendment, "n_covariate_rows": len(rows)}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

