#!/usr/bin/env python3
"""Independently validate the HMSMA stricter proxy-mask sensitivity."""

from __future__ import annotations

import argparse
import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

from spatial_resource_lib import RESOURCE_RELEASE_ID, SpatialResourceError, sha256_file, write_tsv


SENSITIVITY_RELEASE_ID = "hmsma-mask-sensitivity-candidate-2026-08-11-r2"
THRESHOLDS = {500, 1000, 2000}
PROGRAMS = {
    "hotspot_hepatocytes_f05c535ae5bbc0b9",
    "hotspot_hepatocytes_48f39dd4d817a10e",
}


def pdf_pages(path: Path) -> int:
    output = subprocess.run(["pdfinfo", str(path)], check=True, capture_output=True, text=True).stdout
    for line in output.splitlines():
        if line.startswith("Pages:"):
            return int(line.split(":", 1)[1].strip())
    raise SpatialResourceError(f"pdfinfo did not report pages: {path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    project = args.project_root.resolve()
    output = args.output_root.resolve()
    if not (output / "BUILD_COMPLETE").is_file() or (output / "READY").exists():
        raise SpatialResourceError("HMSMA mask candidate is absent, incomplete, or already immutable")
    source_manifest = pd.read_csv(output / "source_manifest.tsv", sep="\t", dtype=str)
    if len(source_manifest) != 42 or source_manifest["source_role"].value_counts().to_dict().get("hmsma_array_h5ad") != 35:
        raise SpatialResourceError("HMSMA source manifest must contain 35 arrays and seven contract sources")
    if source_manifest["relative_path"].duplicated().any():
        raise SpatialResourceError("HMSMA source manifest contains duplicate paths")
    for row in source_manifest.itertuples(index=False):
        path = project / row.relative_path
        if not path.is_file() or path.stat().st_size != int(row.bytes) or sha256_file(path) != row.sha256:
            raise SpatialResourceError(f"HMSMA sensitivity source drift: {row.relative_path}")
    build = pd.read_csv(output / "BUILD_COMPLETE", sep="\t", dtype=str)
    if len(build) != 1 or build.iloc[0]["source_manifest_sha256"] != sha256_file(output / "source_manifest.tsv"):
        raise SpatialResourceError("BUILD_COMPLETE does not seal the HMSMA source manifest")
    detail = pd.read_csv(output / "data/hmsma_proxy_mask_sensitivity.tsv", sep="\t")
    summary = pd.read_csv(output / "data/hmsma_proxy_mask_summary.tsv", sep="\t")
    if len(detail) != 210 or set(detail["min_umi"]) != THRESHOLDS:
        raise SpatialResourceError("mask sensitivity is not the complete 35-by-2-by-3 family")
    if detail["array_id"].nunique() != 35 or set(detail["program_uid"]) != PROGRAMS:
        raise SpatialResourceError("mask sensitivity array or program universe drifted")
    if not np.isfinite(detail["residual_moran_i"]).all():
        raise SpatialResourceError("mask sensitivity contains non-finite Moran statistics")
    if set(detail["inferential_pvalue_authorized"].astype(str).str.upper()) != {"FALSE"}:
        raise SpatialResourceError("mask sensitivity authorizes population inference")
    banned = {"pvalue", "qvalue", "population_pvalue"}.intersection({column.lower() for column in detail.columns})
    if banned:
        raise SpatialResourceError(f"inferential fields leaked into HMSMA sensitivity: {sorted(banned)}")
    pivot = detail.pivot_table(index=["array_id", "program_uid"], columns="min_umi", values="n_spots", aggfunc="first")
    if not ((pivot[500] >= pivot[1000]) & (pivot[1000] >= pivot[2000])).all():
        raise SpatialResourceError("spot counts are not monotone under stricter UMI masks")

    sealed = pd.read_csv(
        project / "Analysis/Multimodal_Program_Projection/candidates" / RESOURCE_RELEASE_ID / "hmsma_label_blind/per_array_program_organization.tsv",
        sep="\t",
    )
    baseline = detail[detail["min_umi"] == 500]
    merged = baseline.merge(sealed, on=["array_id", "program_uid"], suffixes=("_observed", "_sealed"), validate="one_to_one")
    if len(merged) != 70:
        raise SpatialResourceError("baseline reproduction does not contain 70 rows")
    for field in ("n_spots", "n_graph_eligible_spots", "n_directed_edges", "n_genes_measured"):
        if not (merged[f"{field}_observed"].astype(int) == merged[f"{field}_sealed"].astype(int)).all():
            raise SpatialResourceError(f"baseline reproduction drifted in {field}")
    for field in ("retained_l1_weight", "residual_moran_i"):
        if not np.allclose(merged[f"{field}_observed"], merged[f"{field}_sealed"], rtol=0, atol=1e-12):
            raise SpatialResourceError(f"baseline reproduction drifted in {field}")

    expected_rows = []
    baseline_values = baseline[["array_id", "program_uid", "residual_moran_i"]].rename(columns={"residual_moran_i": "baseline"})
    for (uid, threshold), part in detail.groupby(["program_uid", "min_umi"], observed=True):
        joined = part.merge(baseline_values, on=["array_id", "program_uid"], validate="one_to_one")
        values = joined["residual_moran_i"].astype(float)
        base = joined["baseline"].astype(float)
        expected_rows.append({
            "program_uid": uid,
            "min_umi": int(threshold),
            "n_arrays": 35,
            "n_spots_technical": int(part["n_spots"].sum()),
            "median_residual_moran_i": float(values.median()),
            "median_absolute_delta_from_500": float(np.median(np.abs(values - base))),
            "sign_agreement_with_500_fraction": float((np.sign(values) == np.sign(base)).mean()),
            "spearman_rho_with_500": float(values.corr(base, method="spearman")),
        })
    expected = pd.DataFrame(expected_rows)
    checked = summary.merge(expected, on=["program_uid", "min_umi"], suffixes=("_observed", "_expected"), validate="one_to_one")
    if len(checked) != 6:
        raise SpatialResourceError("mask summary does not contain six complete rows")
    for field in ("n_arrays", "n_spots_technical"):
        if not (checked[f"{field}_observed"].astype(int) == checked[f"{field}_expected"].astype(int)).all():
            raise SpatialResourceError(f"mask summary drifted in {field}")
    for field in ("median_residual_moran_i", "median_absolute_delta_from_500", "sign_agreement_with_500_fraction", "spearman_rho_with_500"):
        if not np.allclose(checked[f"{field}_observed"], checked[f"{field}_expected"], rtol=0, atol=1e-12):
            raise SpatialResourceError(f"mask summary drifted in {field}")

    panel = output / "panels/figS_hmsma_proxy_mask_sensitivity.pdf"
    if pdf_pages(panel) != 1:
        raise SpatialResourceError("HMSMA mask panel is not one page")
    if re.search(rb"/Subtype\s*/Type3\b", panel.read_bytes()):
        raise SpatialResourceError("HMSMA mask panel contains Type 3 fonts")
    write_tsv(
        output / "manifest.tsv",
        ("sensitivity_release_id", "spatial_release_id", "artifact", "relative_path", "bytes", "sha256"),
        [
            {
                "sensitivity_release_id": SENSITIVITY_RELEASE_ID,
                "spatial_release_id": RESOURCE_RELEASE_ID,
                "artifact": name,
                "relative_path": path.relative_to(output).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
            for name, path in (
                ("detail", output / "data/hmsma_proxy_mask_sensitivity.tsv"),
                ("summary", output / "data/hmsma_proxy_mask_summary.tsv"),
                ("panel", panel),
            )
        ],
    )
    write_tsv(
        output / "READY",
        ("sensitivity_release_id", "spatial_release_id", "status", "n_rows", "n_panels", "manifest_sha256", "source_manifest_sha256", "population_inference_authorized", "canonical_figure_written"),
        [{
            "sensitivity_release_id": SENSITIVITY_RELEASE_ID,
            "spatial_release_id": RESOURCE_RELEASE_ID,
            "status": "validated_hmsma_proxy_mask_sensitivity_awaiting_adjudication",
            "n_rows": 210,
            "n_panels": 1,
            "manifest_sha256": sha256_file(output / "manifest.tsv"),
            "source_manifest_sha256": sha256_file(output / "source_manifest.tsv"),
            "population_inference_authorized": "FALSE",
            "canonical_figure_written": "FALSE",
        }],
    )
    print("PASS HMSMA 35-array by two-program by three-mask sensitivity")


if __name__ == "__main__":
    main()
