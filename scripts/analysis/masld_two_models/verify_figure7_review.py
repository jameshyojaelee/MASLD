#!/usr/bin/env python3
"""Independently check source identity and saved-result Figure 7 review tables."""
import argparse
import hashlib
import json
import re
import subprocess
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[3]
BENCH = ROOT / "Analysis/MASLD_Model_Benchmark"
SOURCE = BENCH / "executions/encoder-exposure-step4-20260922T175715Z/sources"
EVALUATION = BENCH / "executions/encoder-field-score-21314196/evaluation"
MODELS = ("geneformer", "uce", "scimilarity", "scprint", "scgpt", "transcriptformer")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def library_identity(out):
    manifest = ROOT / "Analysis/SingleCell/continual_integration/candidates/preflight-2026-08-13T224413Z-20206577/contract/cell_manifest.tsv.gz"
    counts = Counter()
    for frame in pd.read_csv(manifest, sep="\t", usecols=["dataset", "library_id", "donor_id"], chunksize=200000):
        part = frame[frame.dataset == "Liver_Atlas"]
        counts.update({key: int(value) for key, value in part.groupby(["library_id", "donor_id"]).size().items()})
    guilliams = pd.read_csv(SOURCE / "sra_runinfo_SRP352824.csv", dtype=str).set_index("Run")
    macparland = pd.read_csv(SOURCE / "sra_runinfo_SRP149989.csv", dtype=str).set_index("Run")
    geo = json.loads((SOURCE / "geo_series_esummary.json").read_text())["result"]
    titles = {sample["accession"]: sample["title"] for uid in geo["uids"] for sample in geo[uid].get("samples", [])}
    public = min(geo[uid]["pdat"].replace("/", "-") for uid in geo["uids"] if geo[uid]["accession"] in {"GSE192740", "GSE192741", "GSE192742"})
    rows = []
    for (library, donor), cells in sorted(counts.items()):
        if library not in guilliams.index or library in macparland.index:
            raise ValueError("Library does not identify Guilliams exclusively: " + library)
        record = guilliams.loc[library]
        title = titles[record.SampleName]
        donor = donor.removeprefix("Liver_Atlas_")
        match = re.search(r"\bH\d{2}\b", title)
        if record.ScientificName != "Homo sapiens" or match is None or match.group() != donor:
            raise ValueError("Source donor/organism disagreement: " + library)
        rows.append({"library_id": library, "resource_donor": donor, "cells": cells,
                     "gsm": record.SampleName, "source_title": title})
    frame = pd.DataFrame(rows)
    frame.to_csv(out / "independent_liver_identity.tsv", sep="\t", index=False)
    catalog = {item["dataset_id"]: item for item in json.loads((SOURCE / "cxg_current_datasets.json").read_text())}
    resource = set(frame.resource_donor)
    guilliams_donors = set(catalog["e84f2780-51e8-4cfa-8aa0-13bbfef677c7"]["donor_id"])
    macparland_donors = set(catalog["a43aa46b-bd16-47fe-bc3e-19a052624e79"]["donor_id"])
    if len(frame) != 42 or len(resource) != 19 or resource != guilliams_donors or resource & macparland_donors:
        raise ValueError("Independently reconstructed donor roster differs")
    return {"libraries": len(frame), "donors": len(resource), "cells": int(frame.cells.sum()),
            "public_date": public, "Guilliams_donor_matches": len(resource & guilliams_donors),
            "MacParland_donor_matches": len(resource & macparland_donors),
            "manifest_sha256": digest(manifest)}


def bootstrap_summary(out):
    table = pd.read_csv(out / "clean_study_paired_gains.tsv", sep="\t")
    points = pd.read_csv(EVALUATION / "study_metrics.tsv", sep="\t").set_index(["model_id", "held_out_study"])
    maximum = 0.0
    with np.load(EVALUATION / "bootstrap_distributions.npz", allow_pickle=False) as saved:
        for row in table.itertuples():
            if row.n_clean_studies == 0:
                if not np.isnan(row.estimate):
                    raise ValueError("Unsupported Geneformer clean-study score")
                continue
            studies = row.clean_studies.split(",")
            model = np.column_stack([saved[row.model_id + "::study_macro_f1::" + study] for study in studies]).mean(axis=1)
            reference = np.column_stack([saved[row.reference_id + "::study_macro_f1::" + study] for study in studies]).mean(axis=1)
            independently_paired = model - reference
            if independently_paired.shape != (10000,):
                raise ValueError("Bootstrap draw count changed")
            estimates = [points.loc[(row.model_id, study), "donor_class_balanced_macro_f1"] for study in studies]
            references = [points.loc[(row.reference_id, study), "donor_class_balanced_macro_f1"] for study in studies]
            expected = [np.mean(estimates) - np.mean(references), *np.quantile(independently_paired, [.025, .5, .975])]
            difference = float(np.max(np.abs(np.array(expected) - [row.estimate, row.lower, row.median, row.upper])))
            maximum = max(maximum, difference)
            if difference > 1e-12:
                raise ValueError("Clean-study paired bootstrap reconstruction differs")
    return {"comparisons": len(table), "maximum_absolute_difference": maximum,
            "supported_numeric_comparisons": int(table.estimate.notna().sum()),
            "draws": 10000, "new_models_fit": False}


def panels(out):
    gains = pd.read_csv(out / "framework7/tables/panel_7f_gain.tsv", sep="\t")
    if not gains.exposure_state.eq("clean_study_subset").all():
        raise ValueError("Clean-study gain retains a global checkpoint exposure label")
    crosses = (gains.lo <= 0) & (gains.hi >= 0)
    if not gains.crosses_zero.eq(crosses).all():
        raise ValueError("Clean-study interval classification differs")
    expected_verdict = np.where(crosses, "interval crosses zero", "interval excludes zero")
    if not gains.verdict.eq(expected_verdict).all():
        raise ValueError("Clean-study interval verdict differs")
    entries = []
    for family, directory in (("framework", out / "framework7/panels"), ("two_model", out / "two_model7")):
        files = sorted(directory.glob("*.pdf"))
        if len(files) != 6:
            raise ValueError("Six panels required per candidate family")
        for path in files:
            pdf = PdfReader(path)
            resources = pdf.pages[0]["/Resources"].get_object()
            fonts = [item.get_object()["/Subtype"] for item in resources.get("/Font", {}).get_object().values()]
            text = subprocess.check_output(["pdftotext", "-layout", str(path), "-"], text=True)
            if len(pdf.pages) != 1 or "/Type3" in fonts:
                raise ValueError("Page or font check failed: " + path.name)
            if family == "two_model":
                original = ROOT / "figures/candidates/figure7_two_models_20260922_21849207" / path.name
                if digest(path) != digest(original):
                    raise ValueError("Unchanged two-model panel differs")
            if family == "framework" and path.name.startswith("fig7D") and "F4 absent" in text:
                raise ValueError("Unresolved GSE193066 primary ceiling retained")
            entries.append({"family": family, "file": path.name, "pages": 1,
                            "Type3_fonts": 0, "sha256": digest(path)})
    return entries


def two_model_numbers(out):
    directory = out / "two_model7"
    manifest = json.loads((directory / "manifest.json").read_text())
    for filename, expected in manifest["source_tables"].items():
        if digest(directory / filename) != expected:
            raise ValueError("Two-model source table hash differs: " + filename)
    sources = {key: Path(value) for key, value in manifest["data_sources"].items()}
    chrom = json.loads((sources["chromatin"] / "results.json").read_text())
    matched = json.loads((sources["within_histology"] / "results.json").read_text())
    secondary = json.loads((sources["secondary_review"] / "results.json").read_text())
    table = pd.read_csv(directory / "Figure7C_source.tsv", sep="\t")
    maximum = 0.0
    for row in table.itertuples():
        if row.model == "interaction_minus_additive":
            sequence = secondary["sequence"]
            expected = [sequence["interaction_minus_additive_original_skill"],
                        *sequence["interaction_minus_additive_original_skill_ci95"]]
        else:
            source = matched if row.model.startswith("matched_") else chrom
            key = row.model.removeprefix("matched_")
            contrast = source["contrasts"][key]
            expected = [contrast["delta_skill"], *contrast["ci95"]]
        difference = float(np.max(np.abs(np.array(expected) - [row.estimate, row.ci_low95, row.ci_high95])))
        maximum = max(maximum, difference)
        if difference > 1e-12:
            raise ValueError("Two-model Figure 7C saved-result comparison differs")
    allele = json.loads((sources["allele"] / "results.json").read_text())
    table = pd.read_csv(directory / "Figure7D_source.tsv", sep="\t")
    for row in table.itertuples():
        if abs(row.four_fold_pooled_MSE - allele["MSE"][row.model]) > 1e-12:
            raise ValueError("Two-model Figure 7D pooled MSE differs")
    mse_rows = len(table)
    people = pd.concat([pd.read_csv(sources[key] / "acquisition_per_participant.tsv.gz", sep="\t")
                        for key in ("acquisition", "histology_sampling")], ignore_index=True)
    people["policy"] = people.policy.replace({"histology": "histology_equal_count_balancing"})
    repeats = people.groupby(["budget_fraction", "policy", "participant_id"])[["base_SSE", "refit_SSE"]].mean()
    reductions = ((repeats.base_SSE - repeats.refit_SSE) / 1000).groupby(level=[0, 1]).mean()
    acquisition = pd.read_csv(directory / "Figure7E_source.tsv", sep="\t")
    for row in acquisition.itertuples():
        if abs(row.participant_pooled_MSE_reduction - reductions.loc[row.budget_fraction, row.policy]) > 1e-12:
            raise ValueError("Two-model Figure 7E participant-pooled reduction differs")
    a = pd.read_csv(directory / "Figure7A_source.tsv", sep="\t")
    b = pd.read_csv(directory / "Figure7B_source.tsv", sep="\t")
    if a.participant_id.nunique() != 99 or a.region_index.nunique() != 24 or len(a) != 99 * 24:
        raise ValueError("Two-model Figure 7A display census differs")
    if len(b) != 1500 or not set(b.fold).issubset({1, 2, 3, 4}):
        raise ValueError("Two-model Figure 7B development display differs")
    return {"source_tables_hash_checked": 6, "panel_C_contrasts": 6,
            "panel_C_maximum_absolute_difference": maximum,
            "panel_D_pooled_MSE_rows": mse_rows, "panel_E_reconstructed_reductions": len(acquisition),
            "panel_A_participants": 99,
            "panel_A_regions": 24, "panel_B_development_variants": 1500,
            "fold0_outcome_rows_parsed": 0}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    out = args.out.resolve()
    report = {"status": "pass", "identity": library_identity(out), "bootstrap": bootstrap_summary(out),
              "panels": panels(out), "two_model_numbers": two_model_numbers(out),
              "fold0_outcome_rows_parsed": 0, "new_models_fit": False}
    for model in MODELS:
        directory = BENCH / "config/artifacts/models" / model
        crosswalk = json.loads((directory / "development_crosswalk.json").read_text())
        audit = json.loads((directory / "exposure_audit.json").read_text())
        if crosswalk["findings"]["Liver_Atlas"]["public_date"] != "2022-01-11":
            raise ValueError("Wrong current Guilliams public date")
        if audit["development_crosswalk_record"]["sha256"] != digest(directory / "development_crosswalk.json"):
            raise ValueError("Current exposure reader digest differs")
    cells = pd.read_csv(out / "encoder_study_exposure.tsv", sep="\t")
    report["exposure_counts"] = dict(Counter(cells.recomputed_disposition))
    if report["exposure_counts"] != {"clean": 54, "confounded_by_exposure": 22, "exposure_unresolved": 8}:
        raise ValueError("Exposure-cell counts differ")
    (out / "verification.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"status": report["status"], "identity": report["identity"], "bootstrap": report["bootstrap"]}))


if __name__ == "__main__":
    main()
