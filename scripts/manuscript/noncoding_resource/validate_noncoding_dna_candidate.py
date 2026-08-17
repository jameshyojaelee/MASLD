#!/usr/bin/env python3
"""Validate a noncoding-DNA candidate and its optional figure package."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--figures", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def read_tsv(path: Path) -> list[dict[str, str]]:
    require(path.is_file(), f"missing table: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_manifest(root: Path) -> None:
    rows = read_tsv(root / "output_manifest.tsv")
    for row in rows:
        path = root / row["relative_path"]
        require(path.is_file(), f"manifest artifact absent: {path}")
        require(path.stat().st_size == int(row["size_bytes"]), f"size drift: {path}")
        digest_field = "sha256" if "sha256" in row else "md5"
        if digest_field == "sha256":
            observed = sha256_file(path)
        else:
            observed = hashlib.md5(path.read_bytes()).hexdigest()  # noqa: S324
        require(observed == row[digest_field], f"checksum drift: {path}")


def validate_candidate(root: Path) -> dict[str, object]:
    require(not (root / "REJECTED.json").exists(), "candidate is explicitly rejected")
    validate_manifest(root)
    members = read_tsv(root / "credible_set_members_annotated.tsv")
    architecture = read_tsv(root / "credible_set_pip_architecture.tsv")
    context = read_tsv(root / "credible_set_regulatory_context.tsv")
    verdict = read_tsv(root / "genetics_noncoding_verdict.tsv")
    execution = json.loads((root / "execution_manifest.json").read_text())
    require(execution["corrected_coloc_used"] is False, "partial COLOC was used")
    require(
        execution["target_gene_inference_included"] is False, "target inference leaked"
    )
    require(len(architecture) == 6467, "credible-set census drift")
    require(len(members) == 42480, "credible-set member census drift")
    passed = [row for row in architecture if row["pip_sum_gate_passed"] == "true"]
    require(len(passed) == 6415, "PIP-sum gate census drift")
    for row in passed:
        masses = [
            float(row[field])
            for field in (
                "protein_altering_pip_mass",
                "canonical_splice_pip_mass",
                "synonymous_or_utr_pip_mass",
                "other_noncoding_pip_mass",
                "unresolved_pip_mass",
            )
        ]
        require(abs(sum(masses) - 1) <= 1e-8, "consequence PIP masses do not sum")
    aggregate_unresolved = sum(
        float(row["unresolved_pip_mass"]) for row in passed
    ) / len(passed)
    require(aggregate_unresolved <= 0.05, "complete consequence annotation gate failed")
    metrics = {row["metric"]: row for row in verdict}
    require(
        metrics["main_figure_architecture_eligible"]["passed"] == "true",
        "main architecture gate failed",
    )
    require(
        metrics["target_transcript_biotype_available"]["passed"] == "false",
        "COLOC target gate unexpectedly open",
    )
    require(
        len(context) == len(architecture), "context/architecture cardinality mismatch"
    )
    joined_fields = set(members[0])
    require(
        "PP.H4.susie" not in joined_fields
        and "coloc" not in " ".join(joined_fields).lower(),
        "COLOC field leaked into DNA table",
    )
    pnpla3 = [row for row in members if row["variant_id"] == "22:44324727:C:G"]
    require(
        pnpla3
        and all(row["consequence_category"] == "protein_altering" for row in pnpla3),
        "PNPLA3 sentinel is not protein-altering",
    )
    tm6sf2 = [row for row in members if row["variant_id"] == "19:19379549:T:C"]
    require(
        tm6sf2
        and all(row["consequence_category"] == "protein_altering" for row in tm6sf2),
        "TM6SF2 sentinel is not protein-altering",
    )
    return {
        "n_credible_sets": len(architecture),
        "n_pip_gate_passed": len(passed),
        "n_members": len(members),
        "aggregate_annotated_pip_mass": 1 - aggregate_unresolved,
        "corrected_coloc_used": False,
    }


def validate_figures(root: Path) -> dict[str, object]:
    validate_manifest(root)
    expected = {
        "fig2_credible_set_pip_architecture.pdf",
        "fig2_noncoding_regulatory_context.pdf",
        "fig4_noncoding_variant_accessibility.pdf",
        "figS_high_pip_variant_lineage_specificity.pdf",
        "figS_credible_set_noncoding_distribution.pdf",
    }
    observed = {path.name for path in root.glob("*.pdf")}
    require(observed == expected, f"figure family drift: {sorted(observed)}")
    for path in sorted(root.glob("*.pdf")):
        info = subprocess.run(
            ["pdfinfo", str(path)], check=True, capture_output=True, text=True
        ).stdout
        require("Pages:           1" in info, f"not single-page: {path}")
        subprocess.run(
            ["gs", "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=nullpage", str(path)],
            check=True,
        )
    arch = read_tsv(root / "fig2_credible_set_pip_architecture_source.tsv")
    require(len(arch) == 10, "architecture source must be 2 scopes x 5 classes")
    for scope in {row["trait_scope"] for row in arch}:
        total = sum(
            float(row["mean_pip_mass"]) for row in arch if row["trait_scope"] == scope
        )
        require(
            abs(total - 1) <= 1e-8, f"architecture mean masses do not sum for {scope}"
        )
    regulatory = read_tsv(root / "fig2_noncoding_regulatory_context_source.tsv")
    require(len(regulatory) == 10, "regulatory context family drift")
    accessibility = read_tsv(root / "fig4_noncoding_variant_accessibility_source.tsv")
    require(
        len(accessibility) == 12, "accessibility family must be 2 scopes x 6 lineages"
    )
    for rows in (regulatory, accessibility):
        require(
            all(0 <= float(row["value"]) <= 1 for row in rows),
            "proportion outside [0,1]",
        )
    specificity = read_tsv(
        root / "figS_high_pip_variant_lineage_specificity_source.tsv"
    )
    require(
        len({row["variant_id"] for row in specificity}) == 30,
        "variant heatmap must show 30 variants",
    )
    require(len(specificity) == 180, "variant heatmap must show 30 x 6 cells")
    return {"pdf_count": len(expected), "source_table_count": len(expected)}


def main() -> None:
    args = parse_args()
    require(not args.output.exists(), f"validation output exists: {args.output}")
    result = {"status": "pass", "candidate": validate_candidate(args.candidate)}
    if args.figures is not None:
        result["figures"] = validate_figures(args.figures)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
