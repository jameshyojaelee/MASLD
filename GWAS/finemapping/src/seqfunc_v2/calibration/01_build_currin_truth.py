#!/usr/bin/env python3
"""Build authoritative Currin caQTL positive/background tables for calibration."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path


def opener(path: Path, mode: str):
    return gzip.open(path, mode, newline="") if path.suffix == ".gz" else path.open(mode, newline="")


def parse_variant(value: str) -> tuple[str, int, str, str]:
    chrom, pos, ref, alt = value.split(":", 3)
    return chrom, int(pos), ref.upper(), alt.upper()


def variant_class(ref: str, alt: str) -> str:
    if len(ref) == len(alt) == 1:
        return "SNV"
    if len(ref) > len(alt):
        return "deletion"
    if len(ref) < len(alt):
        return "insertion"
    return "MNV"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(16 * 1024 * 1024):
            h.update(chunk)
    return h.hexdigest()


def write_gzip_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with gzip.open(tmp, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp, path)


def positive_rows(path: Path, scale: str) -> tuple[list[dict], int]:
    out = []
    excluded = 0
    with opener(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            try:
                chrom, pos, ref, alt = parse_variant(row["lead_variant_ID"])
                beta = float(row["beta"])
            except (KeyError, ValueError):
                excluded += 1
                continue
            ea = row.get("EA", "").upper()
            nea = row.get("NEA", "").upper()
            if scale == "1kb":
                if {ea, nea} != {ref, alt}:
                    excluded += 1
                    continue
                if ea == alt:
                    beta_alt = beta
                    orientation = "effect_alt"
                elif ea == ref:
                    beta_alt = -beta
                    orientation = "effect_ref_flipped_to_alt"
                else:
                    excluded += 1
                    continue
            else:
                # The 1-Mb source omits EA/NEA. Keep source beta but do not use it
                # for signed calibration; magnitude/discrimination remain valid.
                beta_alt = ""
                orientation = "effect_allele_not_released"
            peak_start = int(row["peak_start_hg38"])
            peak_stop = int(row["peak_stop_hg38"])
            out.append(
                {
                    "variant_id_hg38": f"{chrom}:{pos}:{ref}:{alt}",
                    "chr": chrom,
                    "pos_hg38": pos,
                    "ref": ref,
                    "alt": alt,
                    "variant_class": variant_class(ref, alt),
                    "peak_id": row["peak_ID"],
                    "peak_start_hg38": peak_start,
                    "peak_stop_hg38": peak_stop,
                    "peak_center_hg38": (peak_start + peak_stop) // 2,
                    "distance_to_peak_center": pos - ((peak_start + peak_stop) // 2),
                    "positive": 1,
                    "source_scale": scale,
                    "beta_source": beta,
                    "beta_alt": beta_alt,
                    "effect_orientation": orientation,
                    "q_value": row.get("q_val", ""),
                    "p_nominal": row.get("p_val_nominal", ""),
                    "maf": row.get("all_MAF", row.get("lead_MAF", "")),
                    "imputation_r2": row.get("lead_imputation_r2", ""),
                    "block_1mb": f"{chrom}:{pos // 1_000_000}",
                    "tested_status": "source_significant_lead",
                }
            )
    return out, excluded


def negative_rows(path: Path) -> tuple[list[dict], int]:
    out = []
    excluded = 0
    with opener(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            try:
                chrom, pos, ref, alt = parse_variant(row["variant_ID"])
            except (KeyError, ValueError):
                excluded += 1
                continue
            out.append(
                {
                    "variant_id_hg38": f"{chrom}:{pos}:{ref}:{alt}",
                    "chr": chrom,
                    "pos_hg38": pos,
                    "ref": ref,
                    "alt": alt,
                    "variant_class": variant_class(ref, alt),
                    "peak_id": row["peak"],
                    "positive": 0,
                    "source_scale": "non_caPeak_background",
                    "beta_source": "",
                    "beta_alt": "",
                    "effect_orientation": "not_applicable",
                    "q_value": "",
                    "p_nominal": "",
                    "maf": "",
                    "imputation_r2": "",
                    "block_1mb": f"{chrom}:{pos // 1_000_000}",
                    "tested_status": "source_tested_background_beta_adjusted_p_gt_0.5",
                }
            )
    return out, excluded


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(os.environ.get("MASLD_PROJECT_ROOT", Path(__file__).resolve().parents[5])),
    )
    args = parser.parse_args()
    root = args.root.resolve()
    src = root / "GWAS/finemapping/data/seqfunc_external/currin2025_caqtl_v1"
    out = root / "GWAS/finemapping/results/seqfunc/chrombpnet_caqtl/v2_truth"
    out.mkdir(parents=True, exist_ok=True)
    files = {
        "positive_1kb": src / "liver_significant_caQTL_leadVariants_1kb_analysis_with_populationAlleleFrequencies.bed.gz",
        "positive_1mb": src / "liver_significant_caQTL_leadVariants_1mb_analysis.bed.gz",
        "negative": src / "background_variants_overlapping_non-caPeaks.bed.gz",
    }
    if missing := [str(path) for path in files.values() if not path.is_file()]:
        raise SystemExit("Missing Currin files:\n" + "\n".join(missing))
    p1, p1_ex = positive_rows(files["positive_1kb"], "1kb")
    pm, pm_ex = positive_rows(files["positive_1mb"], "1mb")
    neg, neg_ex = negative_rows(files["negative"])
    pfields = list(p1[0])
    nfields = list(neg[0])
    write_gzip_tsv(out / "currin_positive_1kb.tsv.gz", p1, pfields)
    write_gzip_tsv(out / "currin_positive_1mb.tsv.gz", pm, pfields)
    write_gzip_tsv(out / "currin_source_background.tsv.gz", neg, nfields)
    summary = {
        "source": "Currin et al. 2025",
        "zenodo_record": 15025748,
        "release": "v1_2025-03-14",
        "genome_build": "GRCh38",
        "primary_truth": "official 1kb significant lead caQTLs",
        "n_positive_1kb": len(p1),
        "n_positive_1kb_snv": sum(r["variant_class"] == "SNV" for r in p1),
        "n_positive_1mb": len(pm),
        "n_source_background": len(neg),
        "excluded_parse_or_orientation": {"1kb": p1_ex, "1mb": pm_ex, "background": neg_ex},
        "signed_primary_orientation": "beta_alt: positive means alternate allele increases accessibility",
        "one_mb_signed_use": False,
        "input_sha256": {key: sha256(path) for key, path in files.items()},
    }
    tmp = out / "truth_contract.json.tmp"
    tmp.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    os.replace(tmp, out / "truth_contract.json")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
