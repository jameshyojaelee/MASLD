"""Shared contracts for all-117 ambient-RNA recalibration."""

from __future__ import annotations

import hashlib
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


PROJECT_ROOT = Path(
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
PROGRAM_RELEASE = (
    PROJECT_ROOT
    / "Analysis/Multimodal_Program_Projection/candidates"
    / "program-context-v2-stage-corrected-candidate-2026-08-13-v3"
)
PROGRAM_ROOT = PROGRAM_RELEASE / "hotspot"
REGISTRY = PROGRAM_ROOT / "program_registry_v2.tsv"
MEMBERSHIP = PROGRAM_ROOT / "program_membership_v2.tsv"
READY = PROGRAM_ROOT / "READY"
DONOR_METADATA = (
    PROGRAM_RELEASE / "metadata/donor_metadata_stage_corrected_donor.tsv"
)
HOTSPOT_ROOT = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/hotspot_modules"
GLOBAL_ATLAS = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad"
OLD_RAW_EXPORT = (
    PROJECT_ROOT
    / "Analysis/SingleCell/candidates"
    / "igfbp7-ambient-sensitivity-2026-08-13T140836Z"
)
SCORE_SCRIPT = PROJECT_ROOT / "Analysis/SingleCell/scripts/hotspot_modules/501_run_hotspot.py"

LINEAGES = ("cholangiocytes", "fibroblasts", "hepatocytes", "macrophages", "tcells")
LINEAGE_TO_LABEL = {
    "cholangiocytes": "Cholangiocytes",
    "fibroblasts": "Fibroblasts",
    "hepatocytes": "Hepatocytes",
    "macrophages": "Macrophages",
    "tcells": "T cells",
}
LABEL_TO_LINEAGE = {value: key for key, value in LINEAGE_TO_LABEL.items()}
PRIMARY_STAGES = {"Healthy": 0.0, "Steatosis": 1.0, "Steatohepatitis": 2.0}
SEED = 20260815

PAIRING_FILES = {
    "GSE244832": PROJECT_ROOT / "data/GSE244832/metadata/donor_pairing.csv",
    "GSE185477": PROJECT_ROOT / "data/GSE185477/metadata/donor_pairing.csv",
    "GSE202379": PROJECT_ROOT / "data/GSE202379/metadata/donor_pairing.csv",
    "GSE136103": PROJECT_ROOT / "data/GSE136103/metadata/donor_pairing.csv",
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def refuse_existing(path: Path) -> None:
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)


def build_sample_to_donor() -> dict[str, str]:
    """Mirror the frozen donor-collapse map without importing a mutable producer."""
    mapping: dict[str, str] = {}
    for dataset, path in PAIRING_FILES.items():
        table = pd.read_csv(path, dtype=str)
        require({"donor_id", "rna_srrs"} <= set(table), f"pairing schema drift: {path}")
        for row in table.itertuples(index=False):
            donor = f"{dataset}_{row.donor_id}"
            for sample in str(row.rna_srrs).split(";"):
                sample = sample.strip()
                if sample:
                    require(
                        sample not in mapping or mapping[sample] == donor,
                        f"sequencing sample maps to multiple donors: {sample}",
                    )
                    mapping[sample] = donor

    run_file = PROJECT_ROOT / "Liver_Atlas/metadata/SraRunTable.csv"
    sample_file = PROJECT_ROOT / "Liver_Atlas/metadata/GSE192740_sampleInfo_scRNAseq.tsv"
    runs = pd.read_csv(run_file, dtype=str, usecols=["Run", "shortfilename"])
    samples = pd.read_csv(
        sample_file,
        sep="\t",
        dtype=str,
        usecols=["characteristics: shortFileName", "title"],
    ).rename(columns={"characteristics: shortFileName": "shortfilename"})
    liver = runs.merge(samples, on="shortfilename", how="inner", validate="many_to_one")
    liver["donor_short"] = liver["title"].map(
        lambda value: match.group(0) if (match := re.search(r"H[0-9]+", str(value))) else None
    )
    liver = liver.dropna(subset=["donor_short"])
    require(liver["Run"].nunique() == 48, "Liver Atlas run census drift")
    require(liver["donor_short"].nunique() == 19, "Liver Atlas donor census drift")
    for sample, donor_short in liver[["Run", "donor_short"]].itertuples(index=False):
        require(sample not in mapping, f"cross-source sample collision: {sample}")
        mapping[sample] = f"Liver_Atlas_{donor_short}"
    return mapping


def load_donor_metadata() -> pd.DataFrame:
    table = pd.read_csv(DONOR_METADATA, sep="\t", dtype=str)
    required = {
        "sample",
        "dataset",
        "disease_stage_coarse",
        "exclude_stage_analysis",
    }
    require(required <= set(table), "current donor metadata schema drift")
    table = table.rename(columns={"sample": "donor"}).copy()
    table["exclude"] = table["exclude_stage_analysis"].str.lower().isin(
        {"true", "t", "1", "yes"}
    )
    require(not table["donor"].duplicated().any(), "donor metadata is not donor-unique")
    table["stage_ordinal"] = table["disease_stage_coarse"].map(PRIMARY_STAGES)
    return table


def bh_adjust(
    values: np.ndarray | pd.Series, family_size: int | None = None
) -> np.ndarray:
    p = np.asarray(values, dtype=float)
    out = np.full(len(p), np.nan, dtype=float)
    finite = np.isfinite(p)
    if not finite.any():
        return out
    observed = p[finite]
    order = np.argsort(observed)
    ranked = observed[order]
    denominator = len(ranked) if family_size is None else int(family_size)
    require(denominator >= len(ranked), "BH family size is smaller than finite tests")
    q = np.minimum.accumulate(
        (ranked * denominator / np.arange(1, len(ranked) + 1))[::-1]
    )[::-1]
    q = np.minimum(q, 1.0)
    restored = np.empty_like(q)
    restored[order] = q
    out[np.flatnonzero(finite)] = restored
    return out


def linear_fit(y: np.ndarray, x: np.ndarray, contrast: np.ndarray) -> dict[str, float | int]:
    """OLS and HC3 inference for an explicit linear contrast."""
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    contrast = np.asarray(contrast, dtype=float)
    require(len(y) == x.shape[0], "model row mismatch")
    require(x.ndim == 2 and len(contrast) == x.shape[1], "model contrast mismatch")
    require(np.isfinite(y).all() and np.isfinite(x).all(), "non-finite model input")
    rank = int(np.linalg.matrix_rank(x))
    require(rank == x.shape[1], "rank-deficient model")
    n, k = x.shape
    df = n - k
    require(df > 0, "non-positive residual degrees of freedom")
    xtxi = np.linalg.inv(x.T @ x)
    coefficients = xtxi @ x.T @ y
    residuals = y - x @ coefficients
    estimate = float(contrast @ coefficients)
    sigma2 = float(residuals @ residuals) / df
    variance = float(contrast @ (sigma2 * xtxi) @ contrast)
    se = float(np.sqrt(max(variance, 0.0)))
    leverage = np.clip(np.diag(x @ xtxi @ x.T), 0.0, 1.0 - 1e-12)
    omega = (residuals / (1.0 - leverage)) ** 2
    hc3 = xtxi @ (x.T @ (omega[:, None] * x)) @ xtxi
    hc3_variance = float(contrast @ hc3 @ contrast)
    hc3_se = float(np.sqrt(max(hc3_variance, 0.0)))

    def pvalue(value: float, standard_error: float) -> float:
        if standard_error == 0:
            return 0.0 if value != 0 else 1.0
        return float(2.0 * stats.t.sf(abs(value / standard_error), df))

    return {
        "beta": estimate,
        "se": se,
        "pvalue": pvalue(estimate, se),
        "hc3_se": hc3_se,
        "hc3_pvalue": pvalue(estimate, hc3_se),
        "ci_low": estimate - stats.t.ppf(0.975, df) * se,
        "ci_high": estimate + stats.t.ppf(0.975, df) * se,
        "hc3_ci_low": estimate - stats.t.ppf(0.975, df) * hc3_se,
        "hc3_ci_high": estimate + stats.t.ppf(0.975, df) * hc3_se,
        "n": int(n),
        "residual_df": int(df),
        "max_leverage": float(leverage.max()),
    }


def stage_fit(frame: pd.DataFrame, score_column: str) -> dict[str, float | int]:
    data = frame[
        (~frame["exclude"])
        & frame["stage_ordinal"].notna()
        & frame[score_column].notna()
    ].copy()
    require(len(data) >= 20, "fewer than 20 donors in stage model")
    require(data["dataset"].nunique() >= 2, "fewer than two datasets in stage model")
    dummies = pd.get_dummies(data["dataset"], drop_first=True, dtype=float)
    design = np.column_stack(
        [
            np.ones(len(data)),
            data["stage_ordinal"].to_numpy(float),
            dummies.to_numpy(float),
        ]
    )
    result = linear_fit(
        data[score_column].to_numpy(float),
        design,
        np.r_[0.0, 1.0, np.zeros(design.shape[1] - 2)],
    )
    result["n_donors"] = int(data["donor"].nunique())
    result["n_datasets"] = int(data["dataset"].nunique())
    result["n_stage_levels"] = int(data["disease_stage_coarse"].nunique())
    return result


def evidence_state(row: pd.Series) -> str:
    required = ["raw_beta", "corrected_beta", "corrected_qvalue", "delta_qvalue"]
    if any(not np.isfinite(float(row[column])) for column in required):
        return "untestable"
    if not bool(row["primary_selected_frozen"]):
        return "not_applicable_not_frozen_disease_association"
    raw = float(row["raw_beta"])
    corrected = float(row["corrected_beta"])
    corrected_q = float(row["corrected_qvalue"])
    delta_q = float(row["delta_qvalue"])
    if np.sign(raw) != np.sign(corrected):
        if corrected_q < 0.05:
            return "reversal_supported"
        if abs(corrected) < abs(raw) and delta_q < 0.05:
            return "ambient_sensitive_attenuation"
        return "attenuation_indeterminate"
    attenuated = abs(corrected) < abs(raw)
    if corrected_q < 0.05:
        return "retained"
    if attenuated and delta_q < 0.05:
        return "ambient_sensitive_attenuation"
    if attenuated:
        return "attenuation_indeterminate"
    return "attenuation_indeterminate"
