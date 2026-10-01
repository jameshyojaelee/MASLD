#!/usr/bin/env python3
"""GNMT genetics-node rule in 05 and its independent re-derivation in 07.

Synthetic tier-1/2 and master tables only; no project data are read.
Run: python tests/test_gnmt_state.py  (or pytest)
"""

from __future__ import annotations

import importlib.util
import math
import os
import tempfile
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parents[1]
GNMT = "ENSG00000124713"
ABF_STUDY = "2021_34841290_NAFLD_EUR"
UNTESTABLE = "susie_untestable_insufficient_shared_posterior"


def load(name: str, candidate_root: Path):
    os.environ["FIG_CAND_ROOT"] = str(candidate_root)
    spec = importlib.util.spec_from_file_location(name.replace(".py", ""), HERE / name)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_inputs(root: Path, susie: float, master_rows: list[tuple[str, str, str]]) -> dict[str, Path]:
    paths = {name: root / name for name in ("gene_level.csv", "tier12.csv", "master.csv", "tiers.tsv")}
    pd.DataFrame({"gene": ["GNMT", "MAT1A"], "ensembl": [GNMT, "ENSG00000151224"]}).to_csv(
        paths["gene_level.csv"], index=False
    )
    finite = not math.isnan(susie)
    pd.DataFrame([{
        "ensembl": GNMT, "coloc_best_susie_pp4": susie, "coloc_best_abf_pp4": 0.705831502782787,
        "driving_gwas": "UKBB_ALT" if finite else ABF_STUDY,
        "driving_trait": "ALT" if finite else "NAFLD",
        "driving_tier": 2 if finite else 1,
        "any_main": True, "any_supp": True, "tier34_only": False,
        "susie_driving_gwas": "UKBB_ALT" if finite else "",
        "abf_driving_gwas": ABF_STUDY, "abf_driving_trait": "NAFLD", "abf_driving_tier": 1,
    }]).to_csv(paths["tier12.csv"], index=False)
    pd.DataFrame(master_rows, columns=["gwas_name", "ensembl", "method"]).assign(
        gene="x", **{"PP.H4.abf": 0.1}
    ).to_csv(paths["master.csv"], index=False)
    pd.DataFrame({
        "study_name": [ABF_STUDY, "UKBB_ALT", "MVP_Platelet_EUR"],
        "trait": ["NAFLD", "ALT", "Platelet"],
        "tier": [1, 2, 3],
        "tier_label": ["direct_MASLD", "liver_enzyme", "other"],
        "placement": ["main", "main", "supp"],
    }).to_csv(paths["tiers.tsv"], sep="\t", index=False)
    return paths


def run_case(susie: float, master_rows: list[tuple[str, str, str]]) -> tuple[dict, tuple]:
    with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR")) as tmp:
        root = Path(tmp)
        paths = write_inputs(root, susie, master_rows)
        build = load("05_build_figure_sources.py", root / "candidate")
        build.COLOC_TIER12 = str(paths["tier12.csv"])
        build.COLOC_GENE_LEVEL = str(paths["gene_level.csv"])
        build.COLOC_MASTER = str(paths["master.csv"])
        build.GWAS_TIER = paths["tiers.tsv"]
        observed = build.gnmt_genetics()
        validate = load("07_validate_figure_candidate.py", root / "candidate")
        expected = validate.expected_gnmt_genetics(
            paths["tier12.csv"], paths["master.csv"], paths["tiers.tsv"], GNMT
        )
    return observed, expected


def test_untestable_pair_is_not_abf_only() -> None:
    # SuSiE NaN, ABF 0.706, one tier-2 SuSiE pair fails the shared-posterior check.
    observed, expected = run_case(math.nan, [
        ("UKBB_ALT", f"{GNMT}.12", UNTESTABLE),
        (ABF_STUDY, GNMT, "abf_fallback"),
        ("UKBB_ALT", "ENSG00000151224", UNTESTABLE),
    ])
    assert observed["state"] == "genetic_untestable_shared_posterior", observed["state"]
    assert observed["estimate"] is None
    assert "untestable_gwas=UKBB_ALT" in observed["auxiliary"]
    assert expected[0] == "genetic_untestable_shared_posterior", expected


def test_untestable_outside_tier12_is_ignored() -> None:
    # The only failing pair is a tier-3 study, so the ABF-only call stands, with the ABF driver.
    observed, expected = run_case(math.nan, [
        ("MVP_Platelet_EUR", GNMT, UNTESTABLE),
        (ABF_STUDY, GNMT, "abf_fallback"),
    ])
    assert observed["state"] == "abf_only_sensitivity", observed["state"]
    assert f"driver={ABF_STUDY}" in observed["auxiliary"]
    assert expected[0] == "abf_only_sensitivity" and expected[2] == "abf_driving_gwas", expected


def test_supported_susie_precedes_untestable_pair() -> None:
    observed, expected = run_case(0.83, [
        ("UKBB_ALT", GNMT, "susie"),
        (ABF_STUDY, GNMT, UNTESTABLE),
    ])
    assert observed["state"] == "supported_susie_coloc", observed["state"]
    assert "driver=UKBB_ALT" in observed["auxiliary"]
    assert expected[0] == "supported_susie_coloc", expected


if __name__ == "__main__":
    for test in (test_untestable_pair_is_not_abf_only, test_untestable_outside_tier12_is_ignored,
                 test_supported_susie_precedes_untestable_pair):
        test()
        print(f"PASS {test.__name__}")
