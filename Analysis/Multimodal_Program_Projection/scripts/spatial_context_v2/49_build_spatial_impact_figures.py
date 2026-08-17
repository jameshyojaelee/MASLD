#!/usr/bin/env python3
# KEY MESSAGE: positive spatial autocorrelation is not program-specific evidence;
# matched-gene calibration separates IGFBP7 organization from BICC1 background.
"""Build a result-led spatial publication candidate from sealed analyses."""

from __future__ import annotations

import argparse
import os
import shutil
import sys
from pathlib import Path

import matplotlib

matplotlib.use("pdf")
matplotlib.rcParams.update(
    {
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
        "figure.titlesize": 6,
        "axes.linewidth": 0.45,
        "xtick.major.width": 0.45,
        "ytick.major.width": 0.45,
    }
)
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch
import numpy as np
import pandas as pd

from spatial_resource_lib import SpatialResourceError, sha256_file


RELEASE_ID = "spatial-impact-figures-candidate-2026-08-12"
R4_RELEASE_ID = "spatial-publication-figures-candidate-2026-08-11-r4"
PROGRAM_RELEASE_ID = "program-context-v2-candidate-2026-08-07"
CYAN = "#007C91"
CYAN_MID = "#68B7C7"
CYAN_LIGHT = "#CDE7EE"
GRAY = "#9E9E9E"
GRAY_LIGHT = "#E6E6E5"
DARK_GRAY = "#4D4D4D"
INK = "#222222"
WHITE = "#FFFFFF"

PROGRAMS = [
    "hotspot_hepatocytes_f05c535ae5bbc0b9",
    "hotspot_hepatocytes_48f39dd4d817a10e",
]
DATASETS = ["GSE192741", "Vu_et_al_2025"]


def save_pdf(fig: plt.Figure, path: Path) -> None:
    fig.savefig(path, format="pdf", facecolor=WHITE, bbox_inches="tight", dpi=400)
    plt.close(fig)


def write_table(table: pd.DataFrame, path: Path) -> None:
    table.to_csv(path, sep="\t", index=False)


def build_unit_summary(unit: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (program_id, dataset), part in unit.groupby(["program_id", "dataset"], observed=True):
        if dataset == "GSE192741":
            values = part.groupby("aggregation_group", observed=True)["residual_moran_i"].mean()
            unit_label = "donors"
        else:
            values = part.set_index("sample_id")["residual_moran_i"]
            unit_label = "arrays"
        rows.append(
            {
                "program_id": program_id,
                "dataset": dataset,
                "n_reporting_units": int(len(values)),
                "n_positive_units": int((values > 0).sum()),
                "positive_fraction": float((values > 0).mean()),
                "reporting_unit_label": unit_label,
            }
        )
    return pd.DataFrame(rows)


def build_profile(project: Path) -> tuple[pd.DataFrame, dict[str, Path]]:
    candidates = project / "Analysis/Multimodal_Program_Projection/candidates"
    paths = {
        "r4": candidates / R4_RELEASE_ID,
        "registry": candidates / PROGRAM_RELEASE_ID / "hotspot/program_registry_v2.tsv",
        "unit": candidates / "spatial-unit-null-diagnostics-candidate-2026-08-11-r2/data/figS_spatial_unit_heterogeneity.tsv",
        "null": candidates / "spatial-unit-null-diagnostics-candidate-2026-08-11-r2/data/figS_spatial_matched_null_calibration.tsv",
        "effects": candidates / "spatial-resource-candidate-2026-08-11/effects/spatial_program_effects.parquet",
        "full": candidates / "spatial-full-composition-sensitivity-candidate-2026-08-11-r2/full_composition_program_results.tsv",
        "full_ready": candidates / "spatial-full-composition-sensitivity-candidate-2026-08-11-r2/READY",
        "reuse": candidates / "spatial-full-composition-sensitivity-candidate-2026-08-11-r2/matched_set_reuse_audit.tsv",
    }
    required = [paths["r4"] / "READY", paths["registry"], paths["unit"], paths["null"], paths["effects"], paths["full"], paths["full_ready"], paths["reuse"]]
    for path in required:
        if not path.is_file():
            raise SpatialResourceError(f"sealed impact input is absent: {path}")

    registry = pd.read_csv(paths["registry"], sep="\t")
    registry = registry[registry["program_uid"].isin(PROGRAMS)][
        ["program_uid", "module_name", "primary_qvalue", "primary_hc3_qvalue", "stability_median", "robust_display", "membership_sha256", "external_outcomes_read"]
    ].rename(columns={"program_uid": "program_id", "module_name": "program_label"})
    if len(registry) != 2 or not registry["robust_display"].astype(bool).all():
        raise SpatialResourceError("impact family is not the exact two frozen robust-display programs")
    if registry["external_outcomes_read"].astype(str).str.upper().ne("FALSE").any():
        raise SpatialResourceError("frozen registry unexpectedly read external spatial outcomes")

    unit = pd.read_csv(paths["unit"], sep="\t")
    unit_summary = build_unit_summary(unit)
    null = pd.read_csv(paths["null"], sep="\t")
    null = null[
        [
            "program_id", "dataset", "residual_moran_i", "residual_padj", "null_mean",
            "q050", "q500", "q950", "centered_residual_moran_i", "n_null",
            "matched_set_sha256", "source_dependent", "uncertainty_semantics",
        ]
    ].rename(columns={"residual_padj": "primary_qvalue_spatial"})
    effects = pd.read_parquet(paths["effects"])
    effects = effects[effects["program_uid"].isin(PROGRAMS) & effects["dataset_id"].isin(DATASETS)][
        ["program_uid", "dataset_id", "evidence_state", "biological_unit", "n_biological", "technical_unit", "n_technical", "source_dependence"]
    ].rename(columns={"program_uid": "program_id", "dataset_id": "dataset", "evidence_state": "primary_evidence_state"})
    full = pd.read_csv(paths["full"], sep="\t")
    full = full[
        [
            "program_id", "dataset", "residual_moran_i", "residual_null_mean", "residual_padj",
            "robust", "within_source_support", "evidence_state", "residualization",
        ]
    ].rename(
        columns={
            "residual_moran_i": "full_residual_moran_i",
            "residual_null_mean": "full_null_mean",
            "residual_padj": "full_qvalue",
            "robust": "full_robust",
            "within_source_support": "full_within_source_support",
            "evidence_state": "full_evidence_state",
        }
    )
    full["full_centered_residual_moran_i"] = full["full_residual_moran_i"] - full["full_null_mean"]

    profile = null.merge(registry, on="program_id", validate="many_to_one")
    profile = profile.merge(unit_summary, on=["program_id", "dataset"], validate="one_to_one")
    profile = profile.merge(effects, on=["program_id", "dataset"], validate="one_to_one")
    profile = profile.merge(full, on=["program_id", "dataset"], validate="one_to_one")
    profile["dataset_label"] = profile["dataset"].map({"GSE192741": "GSE", "Vu_et_al_2025": "Vu*"})
    profile["program_short"] = profile["program_id"].map({PROGRAMS[0]: "IGFBP7", PROGRAMS[1]: "BICC1"})
    profile["within_source_call"] = np.where(profile["primary_evidence_state"].isin(["supported", "source_dependent"]), "supported", "indeterminate")
    profile["display_state"] = profile["within_source_call"] + np.where(profile["dataset"] == "Vu_et_al_2025", "*", "")
    # Raw sign of the UNCALIBRATED statistic. Spelled out as "I > 0" rather than a
    # bare "+": every row reads positive here, and only the matched-null column
    # separates the programs. A bare plus sign beside an evidence state invites
    # exactly the inference this figure exists to refute.
    profile["unit_sign_label"] = (
        profile["n_positive_units"].astype(str) + "/"
        + profile["n_reporting_units"].astype(str) + " "
        + profile["reporting_unit_label"] + " I > 0"
    )
    profile["claim_boundary"] = np.where(
        profile["program_id"] == PROGRAMS[0],
        "exceeds matched-gene spatial background; composition-linked, not cell-autonomous",
        "positive Moran statistic but does not exceed matched-gene background; indeterminate, not negative",
    )
    order = {(PROGRAMS[0], DATASETS[0]): 0, (PROGRAMS[0], DATASETS[1]): 1, (PROGRAMS[1], DATASETS[0]): 2, (PROGRAMS[1], DATASETS[1]): 3}
    profile["display_order"] = [order[(row.program_id, row.dataset)] for row in profile.itertuples(index=False)]
    profile = profile.sort_values("display_order").reset_index(drop=True)
    return profile, paths


def render_logic(profile: pd.DataFrame, path: Path) -> pd.DataFrame:
    summary = pd.DataFrame(
        [
            {
                "stage": "frozen_input",
                "value": "2 robust_display programs",
                "evidence": "both donor-level disease-associated before spatial analysis",
            },
            {
                "stage": "naive_spatial_readout",
                "value": "4/4 aggregate Moran statistics positive",
                "evidence": "positive spatial autocorrelation in both programs and sources",
            },
            {
                "stage": "matched_null_calibration",
                "value": "expression/detection-matched gene null",
                "evidence": "tests excess program organization rather than generic autocorrelation",
            },
            {
                "stage": "resolved_outcome",
                "value": "IGFBP7 supported; BICC1 indeterminate",
                "evidence": "same frozen starting tier, different physical embodiment",
            },
        ]
    )
    if not (profile["residual_moran_i"] > 0).all():
        raise SpatialResourceError("Figure 4A positive-autocorrelation premise is false")
    fig, ax = plt.subplots(figsize=(4.75, 1.55))
    ax.set_xlim(0, 4.75)
    ax.set_ylim(0, 1.55)
    ax.axis("off")
    boxes = [
        (0.04, 0.48, 0.95, 0.58, "2 frozen robust\nprograms", CYAN_LIGHT),
        (1.22, 0.48, 1.02, 0.58, "both have positive\nMoran's I", CYAN_LIGHT),
        (2.48, 0.48, 1.03, 0.58, "matched-gene\nnull", WHITE),
    ]
    for x, y, width, height, label, color in boxes:
        patch = FancyBboxPatch((x, y), width, height, boxstyle="round,pad=0.02,rounding_size=0.035", facecolor=color, edgecolor=GRAY, linewidth=0.45)
        ax.add_patch(patch)
        ax.text(x + width / 2, y + height / 2, label, ha="center", va="center")
    for start, stop in ((0.99, 1.22), (2.24, 2.48)):
        ax.annotate("", xy=(stop, 0.77), xytext=(start, 0.77), arrowprops={"arrowstyle": "-|>", "lw": 0.5, "color": GRAY})
    ax.annotate("", xy=(3.77, 1.07), xytext=(3.53, 0.83), arrowprops={"arrowstyle": "-|>", "lw": 0.5, "color": GRAY})
    ax.annotate("", xy=(3.77, 0.46), xytext=(3.53, 0.71), arrowprops={"arrowstyle": "-|>", "lw": 0.5, "color": GRAY})
    supported = FancyBboxPatch((3.78, 0.83), 0.91, 0.45, boxstyle="round,pad=0.02,rounding_size=0.035", facecolor=CYAN, edgecolor=GRAY, linewidth=0.45)
    indeterminate = FancyBboxPatch((3.78, 0.24), 0.91, 0.45, boxstyle="round,pad=0.02,rounding_size=0.035", facecolor=GRAY_LIGHT, edgecolor=GRAY, linewidth=0.45)
    ax.add_patch(supported)
    ax.add_patch(indeterminate)
    ax.text(4.235, 1.055, "IGFBP7\nsupported", ha="center", va="center", color=WHITE)
    ax.text(4.235, 0.465, "BICC1\nindeterminate", ha="center", va="center", color=INK)
    ax.text(2.995, 0.27, "Does the program exceed\nexpected autocorrelation?", ha="center", va="top", color=DARK_GRAY)
    save_pdf(fig, path)
    return summary


def render_discrimination(profile: pd.DataFrame, path: Path) -> pd.DataFrame:
    labels = [f"{row.program_short} | {row.dataset_label}" for row in profile.itertuples(index=False)]
    y = np.arange(len(profile))[::-1]
    fig, ax = plt.subplots(figsize=(5.15, 2.30))
    ax.axvline(0, color=GRAY, linewidth=0.45, linestyle="--", zorder=0)
    for yy, row in zip(y, profile.itertuples(index=False), strict=True):
        ax.plot([row.q050, row.q950], [yy, yy], color=GRAY, linewidth=3.0, solid_capstyle="butt", zorder=1)
        ax.scatter(row.q500, yy, s=18, marker="D", color=DARK_GRAY, linewidth=0, zorder=2)
        source_dependent = row.dataset == "Vu_et_al_2025"
        ax.scatter(
            row.residual_moran_i,
            yy,
            s=27,
            marker="^" if source_dependent else "o",
            facecolor=WHITE if source_dependent else CYAN,
            edgecolor=CYAN,
            linewidth=0.8,
            zorder=3,
        )
        ax.text(0.305, yy, row.unit_sign_label, ha="left", va="center", color=DARK_GRAY)
        supported = row.within_source_call == "supported"
        ax.text(0.445, yy, row.display_state, ha="right", va="center", color=CYAN if supported else DARK_GRAY)
    # Column headers make the two text columns readable as a sequence: the raw
    # sign is uniformly positive, the matched-null comparison is what splits them.
    ax.text(0.305, 3.78, "before calibration", ha="left", va="center",
            color=DARK_GRAY, style="italic")
    ax.text(0.445, 3.78, "after calibration", ha="right", va="center",
            color=DARK_GRAY, style="italic")
    ax.set_yticks(y, labels)
    ax.set_xlabel("Residual Moran's I")
    ax.set_xlim(-0.02, 0.46)
    ax.set_ylim(-0.55, 4.15)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0, pad=3)
    ax.tick_params(axis="x", length=2, pad=2)
    ax.legend(
        handles=[
            Line2D([0], [0], color=GRAY, linewidth=3, label="matched genes: 5th–95th"),
            Line2D([0], [0], marker="D", color="none", markerfacecolor=DARK_GRAY, markeredgecolor=DARK_GRAY, markersize=3.7, label="null median"),
            Line2D([0], [0], marker="o", color="none", markerfacecolor=CYAN, markeredgecolor=CYAN, markersize=3.8, label="observed program"),
        ],
        frameon=False,
        ncol=3,
        loc="upper left",
        bbox_to_anchor=(0, 1.01),
        borderaxespad=0,
        handlelength=1.5,
        handletextpad=0.4,
        columnspacing=0.8,
    )
    fig.subplots_adjust(left=0.20, right=0.99, top=0.83, bottom=0.23)
    save_pdf(fig, path)
    return profile


def build(project: Path, output: Path) -> Path:
    candidates = project / "Analysis/Multimodal_Program_Projection/candidates"
    if output.parent.resolve() != candidates.resolve() or output.name != RELEASE_ID:
        raise SpatialResourceError("impact output must be the named isolated candidate")
    if output.exists():
        raise SpatialResourceError(f"refusing to overwrite impact candidate: {output}")
    staging = output.parent / f".{output.name}.incomplete.{os.environ.get('SLURM_JOB_ID', os.getpid())}"
    if staging.exists():
        raise SpatialResourceError(f"impact staging path already exists: {staging}")
    panels = staging / "panels"
    data = staging / "data"
    panels.mkdir(parents=True)
    data.mkdir()

    profile, paths = build_profile(project)
    r4_panels = paths["r4"] / "panels"
    r4_data = paths["r4"] / "data"
    for source in sorted(r4_panels.glob("*.pdf")):
        if source.name in {"fig4a_spatial_firewall.pdf", "fig4f_multimodal_program_summary.pdf"}:
            continue
        shutil.copy2(source, panels / source.name)
    for source in sorted(r4_data.glob("*.tsv")):
        if source.name in {"fig4a_spatial_firewall.tsv", "fig4f_multimodal_program_summary.tsv"}:
            continue
        shutil.copy2(source, data / source.name)

    logic = render_logic(profile, panels / "fig4a_spatial_calibration_logic.pdf")
    discrimination = render_discrimination(profile, panels / "fig4f_matched_null_discrimination.pdf")
    write_table(logic, data / "fig4a_spatial_calibration_logic.tsv")
    write_table(discrimination, data / "fig4f_matched_null_discrimination.tsv")
    write_table(profile, data / "spatial_transportability_profile.tsv")

    # Preserve the complete assay/source audit, but move it out of the main result.
    shutil.copy2(r4_panels / "fig4f_multimodal_program_summary.pdf", panels / "figS_multimodal_program_source_matrix.pdf")
    shutil.copy2(r4_data / "fig4f_multimodal_program_summary.tsv", data / "figS_multimodal_program_source_matrix.tsv")

    claim_contract = pd.DataFrame(
        [
            ("allowed", "positive spatial autocorrelation alone does not establish program-specific organization"),
            ("allowed", "two frozen robust disease programs diverge after matched-gene calibration"),
            ("allowed", "IGFBP7 organization is composition-linked but persists after modeled all-cell adjustment"),
            ("allowed", "BICC1 is spatially indeterminate, not spatially negative"),
            ("prohibited", "BICC1 is not spatial"),
            ("prohibited", "a positive per-unit Moran sign is evidence of program-specific organization"),
            ("prohibited", "IGFBP7 organization is cell autonomous or causal"),
            ("prohibited", "Vu provides donor-level replication"),
            ("prohibited", "the spatial result proves a fibrosis mechanism"),
        ],
        columns=["claim_class", "statement"],
    )
    claim_contract.insert(0, "impact_release_id", RELEASE_ID)
    write_table(claim_contract, staging / "impact_claim_contract.tsv")

    input_paths = [
        (paths["r4"] / "READY", "accepted_r4_publication_family"),
        (paths["registry"], "frozen_program_registry"),
        (paths["unit"], "physical_unit_results"),
        (paths["null"], "matched_null_results"),
        (paths["effects"], "accepted_spatial_evidence_states"),
        (paths["full"], "all_cell_composition_sensitivity"),
        (paths["full_ready"], "all_cell_composition_ready"),
        (paths["reuse"], "matched_set_reuse_audit"),
        (Path(__file__).resolve(), "producer"),
        (Path(__file__).resolve().parent / "50_validate_spatial_impact_figures.py", "validator"),
        (Path(__file__).resolve().parent / "run_spatial_impact_figures.sbatch", "execution_wrapper"),
    ]
    source_rows = []
    for source, role in input_paths:
        if not source.is_file():
            raise SpatialResourceError(f"manifested impact source absent: {source}")
        source_rows.append(
            {
                "impact_release_id": RELEASE_ID,
                "source_role": role,
                "relative_path": source.relative_to(project).as_posix(),
                "bytes": source.stat().st_size,
                "sha256": sha256_file(source),
            }
        )
    write_table(pd.DataFrame(source_rows), staging / "source_manifest.tsv")

    panel_rows = []
    for panel in sorted(panels.glob("*.pdf")):
        panel_rows.append(
            {
                "impact_release_id": RELEASE_ID,
                "panel": panel.stem,
                "relative_path": panel.relative_to(staging).as_posix(),
                "bytes": panel.stat().st_size,
                "sha256": sha256_file(panel),
                "canonical_written": "FALSE",
            }
        )
    write_table(pd.DataFrame(panel_rows), staging / "panel_manifest.tsv")
    write_table(
        pd.DataFrame(
            [
                {
                    "impact_release_id": RELEASE_ID,
                    "status": "rendered_awaiting_independent_validation",
                    "n_panels": len(panel_rows),
                    "n_new_main_panels": 2,
                    "n_retained_r4_panels": len(panel_rows) - 2,
                    "main_source_matrix_moved_to_supplement": "TRUE",
                    "new_inferential_family": "FALSE",
                    "canonical_written": "FALSE",
                }
            ]
        ),
        staging / "BUILD_COMPLETE",
    )
    os.replace(staging, output)
    print(f"rendered result-led spatial impact family: {output}", flush=True)
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, default=None)
    args = parser.parse_args()
    project = args.project_root.resolve()
    output = (args.output_root or project / "Analysis/Multimodal_Program_Projection/candidates" / RELEASE_ID).resolve()
    try:
        build(project, output)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
