#!/usr/bin/env python3
# KEY MESSAGE: HMSMA label-blind organization should not depend on the single
# UMI>=500 proxy tissue mask, but physical arrays remain non-inferential units.
"""Build an isolated HMSMA stricter-mask sensitivity candidate."""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import shutil
import tempfile
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("pdf")
matplotlib.rcParams.update({
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 6,
    "axes.titlesize": 6,
    "axes.labelsize": 6,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 6,
})
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from spatial_resource_lib import (
    EXPECTED_REGISTRY_SHA256,
    PROGRAM_RELEASE_ID,
    RESOURCE_RELEASE_ID,
    SpatialResourceError,
    load_frozen_programs,
    sha256_file,
    write_tsv,
)


SENSITIVITY_RELEASE_ID = "hmsma-mask-sensitivity-candidate-2026-08-11-r2"
THRESHOLDS = (500, 1000, 2000)
CYAN = "#007C91"
CYAN_MID = "#68B7C7"
GRAY = "#9E9E9E"
GRAY_LIGHT = "#E6E6E5"
INK = "#222222"
CONFIRMATORY_UIDS = (
    "hotspot_hepatocytes_f05c535ae5bbc0b9",
    "hotspot_hepatocytes_48f39dd4d817a10e",
)


def load_baseline_module(script_root: Path):
    path = script_root / "22_build_hmsma_label_blind.py"
    spec = importlib.util.spec_from_file_location("sealed_hmsma_baseline", path)
    if spec is None or spec.loader is None:
        raise SpatialResourceError(f"cannot load HMSMA baseline producer: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module, path


def natural_array(value: str) -> tuple[str, int]:
    prefix, _, suffix = value.rpartition("_")
    return prefix, int(suffix) if suffix.isdigit() else 10**9


def write_source_manifest(
    staging: Path,
    project: Path,
    files: list[Path],
    baseline_script: Path,
    sealed_root: Path,
    hotspot: Path,
) -> Path:
    sources = [
        (Path(__file__).resolve(), "sensitivity_producer"),
        (Path(__file__).with_name("38_validate_hmsma_mask_sensitivity.py").resolve(), "independent_validator"),
        (baseline_script.resolve(), "sealed_baseline_producer"),
        ((sealed_root / "READY").resolve(), "sealed_hmsma_ready"),
        ((sealed_root / "per_array_program_organization.tsv").resolve(), "sealed_hmsma_baseline"),
        ((hotspot / "READY").resolve(), "frozen_program_ready"),
        ((project / "Analysis/Spatial/results/preprocessed/HRA007511_starsolo/build_summary.json").resolve(), "hmsma_build_summary"),
    ]
    sources.extend((path.resolve(), "hmsma_array_h5ad") for path in files)
    rows = []
    for path, role in sources:
        if not path.is_file():
            raise SpatialResourceError(f"missing HMSMA sensitivity source {role}: {path}")
        rows.append({
            "sensitivity_release_id": SENSITIVITY_RELEASE_ID,
            "source_role": role,
            "relative_path": path.relative_to(project).as_posix(),
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        })
    manifest = staging / "source_manifest.tsv"
    write_tsv(
        manifest,
        ("sensitivity_release_id", "source_role", "relative_path", "bytes", "sha256"),
        rows,
    )
    return manifest


def exact_baseline_check(detail: pd.DataFrame, sealed_path: Path) -> None:
    sealed = pd.read_csv(sealed_path, sep="\t")
    observed = detail[detail["min_umi"] == 500].copy()
    merged = observed.merge(
        sealed,
        on=["array_id", "program_uid"],
        how="outer",
        suffixes=("_sensitivity", "_sealed"),
        indicator=True,
    )
    if len(merged) != 70 or not (merged["_merge"] == "both").all():
        raise SpatialResourceError("UMI>=500 sensitivity family does not match the sealed 35-by-2 baseline")
    exact_fields = ("program_label", "membership_sha256", "n_spots", "n_graph_eligible_spots", "n_directed_edges", "n_genes_measured")
    for field in exact_fields:
        if not (merged[f"{field}_sensitivity"].astype(str) == merged[f"{field}_sealed"].astype(str)).all():
            raise SpatialResourceError(f"UMI>=500 baseline drift in {field}")
    for field in ("retained_l1_weight", "residual_moran_i"):
        if not np.allclose(
            merged[f"{field}_sensitivity"].astype(float),
            merged[f"{field}_sealed"].astype(float),
            rtol=0,
            atol=1e-12,
        ):
            raise SpatialResourceError(f"UMI>=500 baseline numerical drift in {field}")


def summarize(detail: pd.DataFrame) -> pd.DataFrame:
    baseline = detail[detail["min_umi"] == 500][["array_id", "program_uid", "residual_moran_i"]].rename(
        columns={"residual_moran_i": "baseline_moran_i"}
    )
    rows = []
    for (uid, threshold), part in detail.groupby(["program_uid", "min_umi"], observed=True):
        joined = part.merge(baseline, on=["array_id", "program_uid"], how="left")
        current = joined["residual_moran_i"].astype(float)
        base = joined["baseline_moran_i"].astype(float)
        sign_agree = np.sign(current) == np.sign(base)
        rows.append({
            "sensitivity_release_id": SENSITIVITY_RELEASE_ID,
            "spatial_release_id": RESOURCE_RELEASE_ID,
            "program_uid": uid,
            "program_label": str(part["program_label"].iloc[0]),
            "min_umi": int(threshold),
            "n_arrays": int(part["array_id"].nunique()),
            "n_spots_technical": int(part["n_spots"].sum()),
            "median_residual_moran_i": float(current.median()),
            "min_residual_moran_i": float(current.min()),
            "max_residual_moran_i": float(current.max()),
            "median_absolute_delta_from_500": float(np.median(np.abs(current - base))),
            "sign_agreement_with_500_fraction": float(sign_agree.mean()),
            "spearman_rho_with_500": float(current.corr(base, method="spearman")),
            "population_pvalue_authorized": "FALSE",
        })
    return pd.DataFrame(rows).sort_values(["program_uid", "min_umi"]).reset_index(drop=True)


def render(detail: pd.DataFrame, output: Path) -> Path:
    baseline = detail[detail["min_umi"] == 500][["array_id", "program_uid", "residual_moran_i"]].rename(
        columns={"residual_moran_i": "baseline_moran_i"}
    )
    comparisons = (1000, 2000)
    labels = {
        CONFIRMATORY_UIDS[0]: "Stromal ECM (IGFBP7)",
        CONFIRMATORY_UIDS[1]: "Ductular injury (BICC1)",
    }
    fig, axes = plt.subplots(2, 2, figsize=(4.70, 4.25))
    for row_index, uid in enumerate(CONFIRMATORY_UIDS):
        for column_index, threshold in enumerate(comparisons):
            ax = axes[row_index, column_index]
            part = detail[(detail["program_uid"] == uid) & (detail["min_umi"] == threshold)].merge(
                baseline[baseline["program_uid"] == uid], on=["array_id", "program_uid"], how="left"
            )
            x = part["baseline_moran_i"].to_numpy(float)
            y = part["residual_moran_i"].to_numpy(float)
            lower = float(min(x.min(), y.min()))
            upper = float(max(x.max(), y.max()))
            padding = max(0.015, 0.08 * (upper - lower))
            limits = (lower - padding, upper + padding)
            ax.plot(limits, limits, color=GRAY, linewidth=0.6, linestyle="--", zorder=0)
            ax.scatter(x, y, s=12, color=CYAN if row_index == 0 else CYAN_MID, edgecolor="white", linewidth=0.3)
            rho = pd.Series(x).corr(pd.Series(y), method="spearman")
            sign_agree = float((np.sign(x) == np.sign(y)).mean())
            ax.text(0.03, 0.97, f"rho={rho:.2f}; sign={sign_agree:.2f}", transform=ax.transAxes, ha="left", va="top", color=INK)
            ax.set_xlim(limits)
            ax.set_ylim(limits)
            ax.set_aspect("equal", adjustable="box")
            ax.set_xlabel("UMI >=500 residual Moran's I")
            ax.set_ylabel(f"UMI >={threshold} residual Moran's I")
            if row_index == 0:
                ax.set_title(f"stricter proxy mask: {threshold:,} UMIs")
            if column_index == 0:
                ax.text(-0.30, 0.5, labels[uid], transform=ax.transAxes, rotation=90, ha="center", va="center")
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_linewidth(0.35)
            ax.spines["bottom"].set_linewidth(0.35)
            ax.tick_params(width=0.35, length=2)
    fig.text(0.5, 0.005, "35 physical arrays; proxy-mask sensitivity only; donor identity unresolved; no population P values.", ha="center", va="bottom", color="#4D4D4D")
    fig.tight_layout(rect=(0.06, 0.04, 1, 1), h_pad=1.0, w_pad=0.8)
    path = output / "figS_hmsma_proxy_mask_sensitivity.pdf"
    fig.savefig(path, bbox_inches="tight", dpi=400)
    plt.close(fig)
    return path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    project = args.project_root.resolve()
    output = args.output_root.resolve()
    if output.exists():
        raise SpatialResourceError(f"immutable HMSMA mask candidate already exists: {output}")
    source_root = project / "Analysis/Spatial/results/preprocessed/HRA007511_starsolo"
    sealed_root = project / "Analysis/Multimodal_Program_Projection/candidates" / RESOURCE_RELEASE_ID / "hmsma_label_blind"
    if not (sealed_root / "READY").is_file():
        raise SpatialResourceError("sealed HMSMA label-blind baseline is missing")
    summary = json.loads((source_root / "build_summary.json").read_text(encoding="utf-8"))
    if summary.get("min_umi") != 500 or summary.get("n_samples_built") != 35 or summary.get("phenotype_key_available") is not False:
        raise SpatialResourceError(f"HMSMA source contract drift: {summary}")
    files = sorted(source_root.glob("HRA_*.h5ad"), key=lambda path: natural_array(path.stem))
    if len(files) != 35:
        raise SpatialResourceError(f"expected 35 HMSMA arrays, found {len(files)}")
    script_root = Path(__file__).resolve().parent
    baseline_module, baseline_script = load_baseline_module(script_root)
    hotspot = project / "Analysis/Multimodal_Program_Projection/candidates" / PROGRAM_RELEASE_ID / "hotspot"
    programs, weights = load_frozen_programs(hotspot)
    if not set(CONFIRMATORY_UIDS).issubset(programs):
        raise SpatialResourceError("frozen confirmatory program family is incomplete")

    staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.", dir=output.parent))
    try:
        source_manifest = write_source_manifest(
            staging, project, files, baseline_script, sealed_root, hotspot
        )
        rows = []
        for path in files:
            adata = ad.read_h5ad(path)
            try:
                total_all = pd.to_numeric(adata.obs["total_counts"], errors="coerce").to_numpy(float)
                detected_all = pd.to_numeric(adata.obs["n_genes"], errors="coerce").to_numpy(float)
                symbols = baseline_module.gene_symbols(adata)
                for threshold in THRESHOLDS:
                    keep = np.isfinite(total_all) & (total_all >= threshold)
                    if int(keep.sum()) < 50:
                        raise SpatialResourceError(f"{path.stem} has fewer than 50 spots at UMI>={threshold}")
                    matrix = adata.X[keep]
                    total = total_all[keep]
                    detected = detected_all[keep]
                    obs = adata.obs.loc[keep]
                    for uid in CONFIRMATORY_UIDS:
                        residual, measured, retained = baseline_module.score_program(
                            matrix, total, detected, symbols, weights[uid]
                        )
                        statistic, n_graph, n_edges = baseline_module.moran_hex(
                            residual,
                            obs["array_row"].to_numpy(),
                            obs["array_col"].to_numpy(),
                        )
                        rows.append({
                            "sensitivity_release_id": SENSITIVITY_RELEASE_ID,
                            "spatial_release_id": RESOURCE_RELEASE_ID,
                            "program_release_id": PROGRAM_RELEASE_ID,
                            "registry_sha256": EXPECTED_REGISTRY_SHA256,
                            "dataset_id": "HRA007511_HMSMA",
                            "array_id": path.stem,
                            "program_uid": uid,
                            "program_label": programs[uid]["module_name"],
                            "membership_sha256": programs[uid]["membership_sha256"],
                            "min_umi": threshold,
                            "mask_role": "sealed_baseline" if threshold == 500 else "stricter_proxy_mask_sensitivity",
                            "n_spots": int(keep.sum()),
                            "n_graph_eligible_spots": n_graph,
                            "n_directed_edges": n_edges,
                            "n_genes_measured": len(measured),
                            "retained_l1_weight": retained,
                            "residual_moran_i": statistic,
                            "score_unit": "frozen_positive_weight_log1p_CPT_residual_z",
                            "adjustment": "within_array_log_library_size_and_detected_gene_count",
                            "inferential_pvalue_authorized": "FALSE",
                        })
            finally:
                del adata
        detail = pd.DataFrame(rows).sort_values(["array_id", "program_uid", "min_umi"]).reset_index(drop=True)
        exact_baseline_check(detail, sealed_root / "per_array_program_organization.tsv")
        summary_table = summarize(detail)
        data_dir = staging / "data"
        panel_dir = staging / "panels"
        data_dir.mkdir(parents=True)
        panel_dir.mkdir(parents=True)
        detail.to_csv(data_dir / "hmsma_proxy_mask_sensitivity.tsv", sep="\t", index=False)
        summary_table.to_csv(data_dir / "hmsma_proxy_mask_summary.tsv", sep="\t", index=False)
        panel = render(detail, panel_dir)
        write_tsv(
            staging / "BUILD_COMPLETE",
            ("sensitivity_release_id", "spatial_release_id", "status", "thresholds", "n_arrays", "n_programs", "n_rows", "source_manifest_sha256", "baseline_producer_sha256", "baseline_ready_sha256", "panel_sha256", "population_inference_authorized", "canonical_figure_written"),
            [{
                "sensitivity_release_id": SENSITIVITY_RELEASE_ID,
                "spatial_release_id": RESOURCE_RELEASE_ID,
                "status": "rendered_candidate_awaiting_independent_validation",
                "thresholds": ";".join(map(str, THRESHOLDS)),
                "n_arrays": 35,
                "n_programs": 2,
                "n_rows": len(detail),
                "source_manifest_sha256": sha256_file(source_manifest),
                "baseline_producer_sha256": sha256_file(baseline_script),
                "baseline_ready_sha256": sha256_file(sealed_root / "READY"),
                "panel_sha256": sha256_file(panel),
                "population_inference_authorized": "FALSE",
                "canonical_figure_written": "FALSE",
            }],
        )
        os.replace(staging, output)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    print(json.dumps({"release": SENSITIVITY_RELEASE_ID, "rows": len(detail), "thresholds": THRESHOLDS}, indent=2))


if __name__ == "__main__":
    main()
