#!/usr/bin/env python3
"""Prepare two Figure 7 review sets from saved results, without model fitting."""
import argparse
import csv
import hashlib
import json
import os
import re
import shutil
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
BENCH = ROOT / "Analysis/MASLD_Model_Benchmark"
EXPOSURE = BENCH / "executions/encoder-exposure-step4-20260922T175715Z"
FRAMEWORK = BENCH / "executions/fig7_panels_rebuild-20260901T012245Z"
EVALUATION = BENCH / "executions/encoder-field-score-21314196/evaluation"
TWO_MODEL = ROOT / "figures/candidates/figure7_two_models_20260922_21849207"
MODELS = ("geneformer", "uce", "scimilarity", "scprint", "scgpt", "transcriptformer")
UUID = "e84f2780-51e8-4cfa-8aa0-13bbfef677c7"
ALIASES = ["Liver_Atlas", "GSE192742", "GSE192740", "GSE192741",
           "DOI:10.1016/j.cell.2021.12.018", "CELLXGENE:" + UUID]
CLEAN = {"clean_declared", "target_label_unexposed"}
CONFOUNDED = {"encoder_seen", "continual_seen", "reference_only"}


def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1048576), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path, value):
    with path.open("x") as handle:
        json.dump(value, handle, indent=2)
        handle.write("\n")


def replace_object(text, key, value):
    """Replace one JSON member while retaining all unrelated file formatting."""
    match = re.search(r'(?m)^([ \t]*)"' + re.escape(key) + r'"\s*:\s*', text)
    if match is None:
        raise ValueError("JSON member not found: " + key)
    start = match.end()
    _old, end = json.JSONDecoder().raw_decode(text[start:])
    encoded = json.dumps(value, indent=2, ensure_ascii=False)
    encoded = encoded.replace("\n", "\n" + match.group(1))
    return text[:start] + encoded + text[start + end:]


def verify_evidence():
    reconstruction = json.loads((EXPOSURE / "out/encoder_study_exposure_crosswalk.json").read_text())
    for filename, record in reconstruction["sources"].items():
        if sha(EXPOSURE / "sources" / filename) != record["sha256"]:
            raise ValueError("Exposure primary source hash changed: " + filename)
    identity = pd.read_csv(EXPOSURE / "out/liver_atlas_identity_check.tsv", sep="\t")
    if (len(identity) != 42 or identity.library_id.nunique() != 42
            or identity.resource_donor.nunique() != 19
            or not identity.in_SRP352824_GSE192740.all()
            or identity.in_SRP149989_GSE115469.any()
            or not identity.donor_agrees.all()
            or set(identity.organism) != {"Homo sapiens"}):
        raise ValueError("Library identity reconstruction differs")
    return reconstruction


def correct_crosswalks(out, reconstruction, apply_current):
    directory = BENCH / "config/artifacts/models"
    snapshot = out / "config_before"
    snapshot.mkdir()
    proposed = out / "corrected_crosswalks"
    proposed.mkdir()
    changed, digest_map = [], {}
    blocks = reconstruction["checkpoints"]
    for model in MODELS:
        source = directory / model / "development_crosswalk.json"
        old_text = source.read_text()
        old = json.loads(old_text)
        current = dict(old["findings"]["Liver_Atlas"])
        current.update(aliases=ALIASES, public_date="2022-01-11",
                       article_title="Spatial proteogenomics reveals distinct and evolutionarily conserved hepatic macrophage niches")
        if model == "geneformer":
            finding = blocks["geneformer_v2_316m"]["findings"]["Liver_Atlas"]
            current["corpus_findings"] = {
                "Genecorpus-104M": {"declared_match": True, "exposure_state": "encoder_seen",
                    "evidence": finding["evidence"], "matched_records": [
                        {"cells_passed_filtering": record["cells"],
                         "source_key": record["source"] + ":row" + str(record["row"])}
                        for record in finding["matched_records"]]},
                "Genecorpus-30M": {"declared_match": False, "exposure_state": "clean_declared",
                    "matched_records": [],
                    "evidence": "Guilliams GEO public date 2022-01-11 is after the documented June 2021 V1 training cutoff; MacParland corpus rows are a different study."}}
        elif model in {"uce", "scimilarity", "transcriptformer"}:
            key = {"uce": "uce", "scimilarity": "scimilarity_v1_1",
                   "transcriptformer": "transcriptformer_tf_sapiens"}[model]
            finding = blocks[key]["findings"]["Liver_Atlas"]
            current.update({k: finding[k] for k in ("evidence", "exposure_state", "matched_records")})
        elif model == "scprint":
            census = pd.read_csv(EXPOSURE / "sources/census_2023-12-15_datasets.tsv", sep="\t", dtype=str)
            hit = census[census.dataset_id == UUID]
            if len(hit) != 1:
                raise ValueError("Guilliams UUID not unique in frozen scPrint-era Census")
            current.pop("historical_archive_records", None)
            current.update(exposure_state="unknown", evidence="Guilliams dataset " + UUID + " is present in the frozen Census 2023-12-15 inventory; the checkpoint's post-QC training/holdout membership is unpublished and remains unknown.",
                           historical_archive_records=[{"dataset_id": UUID, "census_build": "2023-12-15"}])
        else:
            current["evidence"] = "Guilliams GSE192740/GSE192742 was public before the checkpoint cutoff. Membership in the retired 2023-05-08 weekly Census build remains unavailable; no later build substitutes for it."
        new_text = replace_object(old_text, "Liver_Atlas", current)
        if model == "geneformer":
            interpretation = dict(old["interpretation"])
            interpretation["geneformer_v1_10m"] = "GSE136103 is encoder_seen. Guilliams Liver_Atlas first became public after the June 2021 V1 cutoff and is clean_declared; MacParland rows do not identify this Resource study. Other checkpoint-specific findings and protected GSE289173 labels are unchanged."
            new_text = replace_object(new_text, "interpretation", interpretation)
        if model == "transcriptformer":
            findings = json.loads(new_text)["findings"]
            for study, finding in blocks["transcriptformer_tf_sapiens"]["findings"].items():
                if study == "Liver_Atlas" or study not in findings:
                    continue
                replacement = dict(findings[study])
                replacement.update({k: finding[k] for k in ("evidence", "exposure_state", "matched_records")})
                new_text = replace_object(new_text, study, replacement)
            corpus = dict(old["corpus"])
            declared = reconstruction["sources"]["transcriptformer_pearce2025_ED_table1.xlsx"]
            corpus.update(source_record_sha256=declared["sha256"], source_record_status="declared_644_human_UUID_inventory_resolved_for_seven_Resource_studies_2026-09-22",
                          declared_cell_total=57336190)
            new_text = replace_object(new_text, "corpus", corpus)
            scope = dict(old["crosswalk_scope"])
            scope["scope_limit"] = "Cell-state task only. The seven Resource studies and GSE296875 use the resolved declared UUID inventory. Other pre-archive sources, including protected GSE289173, retain their existing unresolved states."
            new_text = replace_object(new_text, "crosswalk_scope", scope)
            for study in ("GSE256398", "GSE289173"):
                finding = dict(json.loads(new_text)["findings"][study])
                finding["evidence"] = "Not adjudicated by the seven-Resource-study UUID reconstruction. Its prior unknown exposure state is retained; available corpus metadata alone is not an exposure decision for this source."
                new_text = replace_object(new_text, study, finding)
        if model == "uce":
            corpus = dict(old["corpus"])
            corpus["declared_dataset_rows"] = blocks["uce"]["summary"]["rows"]
            new_text = replace_object(new_text, "corpus", corpus)
        new_text = replace_object(new_text, "verified_date", "2026-10-05")
        shutil.copyfile(source, snapshot / (model + "_development_crosswalk.json"))
        destination = proposed / (model + "_development_crosswalk.json")
        destination.write_text(new_text)
        digest_map[sha(source)] = sha(destination)
        if apply_current:
            source.write_text(new_text)
        changed.append(str(source.relative_to(ROOT)))
    for model in ("geneformer", "transcriptformer"):
        source = directory / model / "exposure_audit.json"
        old_text = source.read_text()
        if model == "geneformer":
            interpretation = json.loads(old_text)["interpretation"].replace(
                "The cell-state crosswalk detects declared encoder exposure for GSE136103 and Liver_Atlas in V1 and V2, and GSE185477 in V2.",
                "The corrected cell-state crosswalk detects GSE136103 exposure in V1 and V2, and Guilliams Liver_Atlas plus GSE185477 exposure in V2. Guilliams was first public after the V1 cutoff and is clean_declared there.")
            new_text = replace_object(old_text, "interpretation", interpretation)
        else:
            new_text = replace_object(old_text, "interpretation", "The declared 644-human-UUID inventory is now available and resolves the seven Resource studies: five are clean and GSE185477 plus Guilliams Liver_Atlas are encoder_seen. This reconstruction does not adjudicate protected GSE289173, whose separate prior unknown exposure finding and eligibility remain unchanged. The mutable catalog is used only for identity, not as a replacement corpus.")
            new_text = replace_object(new_text, "remaining_blockers", ["Protected GSE289173 was not adjudicated in the seven-study reconstruction; its prior unknown finding is retained.", "Do not transfer this audit to TF-Exemplar, TF-Metazoa, later releases, or independently tuned checkpoints."])
        new_text = replace_object(new_text, "verified_date", "2026-10-05")
        shutil.copyfile(source, snapshot / (model + "_exposure_audit.json"))
        destination = proposed / (model + "_exposure_audit.json")
        destination.write_text(new_text)
        digest_map[sha(source)] = sha(destination)
        if apply_current:
            source.write_text(new_text)
        changed.append(str(source.relative_to(ROOT)))
    if apply_current:
        # Only current config references are refreshed; historical executions and releases remain intact.
        for _ in range(8):
            new_digests = {}
            for path in sorted((BENCH / "config").rglob("*.json")):
                text = path.read_text()
                refreshed = text
                for previous, current in digest_map.items():
                    refreshed = refreshed.replace('"' + previous + '"', '"' + current + '"')
                if refreshed != text:
                    previous = sha(path)
                    safe_name = str(path.relative_to(BENCH / "config")).replace("/", "__")
                    snapshot_path = snapshot / safe_name
                    if not snapshot_path.exists():
                        snapshot_path.write_text(text)
                    path.write_text(refreshed)
                    new_digests[previous] = sha(path)
                    changed.append(str(path.relative_to(ROOT)))
            if not new_digests:
                break
            digest_map.update(new_digests)
        else:
            raise ValueError("Dependent config digest updates did not settle")
    return sorted(set(changed))


def exposure_summary(out, reconstruction):
    cells = pd.read_csv(EXPOSURE / "out/fig7e_cell_dispositions.tsv", sep="\t")
    counts = Counter(cells.recomputed_disposition)
    if counts != {"clean": 54, "confounded_by_exposure": 22, "exposure_unresolved": 8}:
        raise ValueError("84-cell exposure summary differs")
    cells.to_csv(out / "encoder_study_exposure.tsv", sep="\t", index=False)
    metrics = pd.read_csv(EVALUATION / "study_metrics.tsv", sep="\t")
    points = metrics.set_index(["model_id", "held_out_study"]).donor_class_balanced_macro_f1
    names = {"Geneformer": "geneformer_v2_316m", "TranscriptFormer": "transcriptformer_tf_sapiens",
             "SCimilarity": "scimilarity_v1_1", "UCE-4L": "uce_4l", "UCE-33L": "uce_33l", "HVG+PCA": "hvg_pca_task_native"}
    head_names = {"linear head": "linear", "two-layer MLP head": "two_layer_mlp"}
    rows = []
    with np.load(EVALUATION / "bootstrap_distributions.npz", allow_pickle=False) as draws:
        for (encoder, head), part in cells.groupby(["encoder", "head"], sort=True):
            studies = sorted(part.loc[part.recomputed_disposition == "clean", "held_out_study"])
            model = names[encoder] + "::" + head_names[head]
            for reference in ("hvg_pca_task_native", "scvi_liver_latent"):
                ref = reference + "::" + head_names[head]
                row = {"encoder": encoder, "head": head, "model_id": model, "reference_id": ref,
                       "clean_studies": ",".join(studies), "n_clean_studies": len(studies),
                       "n_clean_donor_records": int(part.loc[part.recomputed_disposition == "clean", "donors"].sum()),
                       "uncertainty": "saved paired 10000-draw study-stratified donor-cluster bootstrap; equal study weighting; conditional on fixed models"}
                if not studies:
                    row.update(estimate=np.nan, lower=np.nan, median=np.nan, upper=np.nan,
                               state="uninterpretable_no_clean_study")
                else:
                    values = np.mean([draws[model + "::study_macro_f1::" + study] - draws[ref + "::study_macro_f1::" + study] for study in studies], axis=0)
                    row.update(estimate=float(np.mean([points[model, study] - points[ref, study] for study in studies])),
                               lower=float(np.quantile(values, .025)), median=float(np.quantile(values, .5)),
                               upper=float(np.quantile(values, .975)), state="clean_study_comparison")
                rows.append(row)
    pd.DataFrame(rows).to_csv(out / "clean_study_paired_gains.tsv", sep="\t", index=False)
    return cells, pd.DataFrame(rows), dict(counts)


def framework_candidate(out, cells, gains):
    out.mkdir()
    shutil.copytree(FRAMEWORK / "tables", out / "tables")
    shutil.copytree(FRAMEWORK / "scripts", out / "scripts")
    (out / "panels").mkdir()
    table = out / "tables/panel_7f_exposure_grid.tsv"
    frame = pd.read_csv(table, sep="\t")
    updated = cells.set_index(["encoder", "head", "held_out_study"])
    short = {"geneformer_v2_316m": "Geneformer", "transcriptformer_tf_sapiens": "TranscriptFormer",
             "scimilarity_v1_1": "SCimilarity", "uce_4l": "UCE-4L", "uce_33l": "UCE-33L", "hvg_pca_task_native": "HVG+PCA"}
    for index, entry in frame.iterrows():
        row = updated.loc[short[entry.block_id], entry.head_label, entry.held_out_study]
        frame.loc[index, ["study_exposure_state", "disposition"]] = [row.recomputed_state, row.recomputed_disposition]
    frame.to_csv(table, sep="\t", index=False)
    table = out / "tables/panel_7f_gain.tsv"
    previous = pd.read_csv(table, sep="\t")
    previous.to_csv(out / "tables/conditional_all_study_panel_7f_gain.tsv", sep="\t", index=False)
    selected = previous[previous.block.str.startswith("b.")].copy()
    clean_rows = gains[gains.reference_id.str.startswith("hvg_pca_task_native")]
    for index, entry in selected.iterrows():
        hit = clean_rows[clean_rows.model_id == entry.model_id]
        if entry.block_id == "scvi_liver_latent":
            # Reverse the same-study HVG-minus-liver contrast, not a pooled donor estimand.
            hit = gains[(gains.model_id == "hvg_pca_task_native::" + entry.head_id)
                        & (gains.reference_id == entry.model_id)].iloc[0]
            selected.loc[index, ["gain", "lo", "hi"]] = [-hit.estimate, -hit.upper, -hit.lower]
            selected.loc[index, "note"] = "Equal-study mean over all seven clean studies; paired saved bootstrap."
            selected.loc[index, "unit_and_n"] = str(hit.n_clean_donor_records) + " donor records, 7 clean studies; equal study weights"
            selected.loc[index, "source"] = str((out.parent / "clean_study_paired_gains.tsv").relative_to(ROOT))
            continue
        if len(hit) != 1:
            raise ValueError("Clean gain row not unique")
        hit = hit.iloc[0]
        selected.loc[index, ["gain", "lo", "hi"]] = [hit.estimate, hit.lower, hit.upper]
        selected.loc[index, "clean_studies"] = hit.n_clean_studies
        relevant = frame[frame.model_id == entry.model_id]
        selected.loc[index, "confounded_studies"] = sum(relevant.disposition == "confounded_by_exposure")
        selected.loc[index, "unresolved_studies"] = sum(relevant.disposition == "exposure_unresolved")
        selected.loc[index, "eligibility"] = hit.state
        selected.loc[index, "unit_and_n"] = str(hit.n_clean_donor_records) + " donor records, " + str(hit.n_clean_studies) + " clean studies; equal study weights"
        selected.loc[index, "note"] = "Equal-study mean over clean studies only; comparator uses the identical study subset. Conditional on saved fitted models."
        selected.loc[index, "source"] = str((out.parent / "clean_study_paired_gains.tsv").relative_to(ROOT))
    selected = selected[selected.gain.notna()].copy()
    selected["crosses_zero"] = (selected.lo <= 0) & (selected.hi >= 0)
    selected["verdict"] = np.where(selected.crosses_zero, "interval crosses zero", "interval excludes zero")
    selected.to_csv(table, sep="\t", index=False)
    model_metrics = pd.read_csv(EVALUATION / "model_metrics.tsv", sep="\t")
    for index, entry in model_metrics.iterrows():
        relevant = frame[frame.model_id == entry.model_id]
        if relevant.empty:
            continue
        per = relevant.groupby("outer_fold").disposition.apply(lambda x: "confounded_by_exposure" if "confounded_by_exposure" in set(x) else "exposure_unresolved" if "exposure_unresolved" in set(x) else "clean")
        model_metrics.loc[index, "clean_folds"] = sum(per == "clean")
    model_metrics.to_csv(out / "tables/model_metrics.tsv", sep="\t", index=False)
    script = out / "scripts/fig7f_liver_encoder_against_the_field.R"
    text = script.read_text().replace('MM_REL <- "Analysis/MASLD_Model_Benchmark/executions/encoder-field-score-21314196/evaluation/model_metrics.tsv"', 'MM_REL <- file.path(Sys.getenv("FIG7_EXEC_DIR"), "tables", "model_metrics.tsv")')
    text = text.replace('file.path(ROOT, MM_REL)', 'MM_REL')
    text = text.replace('". The two best-pooled encoders, ", lev_all[1], " and ", lev_all[2],\n        ", have ZERO clean held-out studies, so their lead cannot be read as skill. ",', '". Pooled scores mix exposure states and are descriptive. TranscriptFormer now has five clean held-out studies; Geneformer has none. ",')
    text = text.replace('"gain over HVG + PCA"', '"clean-study gain over HVG + PCA"')
    text = text.replace('scale_x_continuous(breaks = c(0, 0.02, 0.04, 0.06))', 'scale_x_continuous(n.breaks = 4)')
    text = text.replace('"Bottom: paired gain in donor-class-balanced macro F1 over HVG + 50-d PCA of raw ",', '"Bottom: paired gain on each encoder\'s clean studies only; the HVG + 50-d PCA comparator uses exactly the same studies. Equal-study means of donor-class-balanced macro F1 from saved predictions. Raw ",')
    text = text.replace('"; with a two-layer MLP head the same encoders gain only "', '"; with a two-layer MLP head the corresponding gains range "')
    text = text.replace('"counts carrying the identical head; point = estimate, whisker = 95% donor cluster "', '"counts carrying the identical head; point = estimate, whisker = nominal 95% donor cluster "')
    text = text.replace('"bootstrap interval; magenta the interval excludes zero, gray it crosses zero. With "', '"bootstrap interval; magenta is above zero, blue below zero, and gray crosses zero. With "')
    script.write_text(text)
    bounds = pd.read_csv(out / "tables/panel_7d_bounds.tsv", sep="\t")
    withheld = bounds[bounds.row.str.contains("GSE193066", na=False)].copy()
    withheld["evidence_state"] = "conditional_identity_unresolved_not_primary"
    withheld.to_csv(out / "tables/conditional_GSE193066_bounds.tsv", sep="\t", index=False)
    bounds = bounds[~bounds.row.str.contains("GSE193066", na=False)]
    bounds.to_csv(out / "tables/panel_7d_bounds.tsv", sep="\t", index=False)
    script = out / "scripts/fig7d_measured_bounds_and_join_structure.R"
    text = script.read_text()
    text = re.sub(r'(?m)^\s*"GSE193066 fibrosis stage 0-3, cross-sectional, n=106 participants"[^\n]*\n', '', text)
    script.write_text(text)
    (out / "REVIEW_NOTES.md").write_text("Candidate only. Current working Figure 7 and all model releases are unchanged.\n\nFigure 7F uses corrected Guilliams identity and the resolved TranscriptFormer UUID inventory: 54 clean, 22 confounded, 8 unresolved encoder/head/study cells. Lower comparisons use each encoder's clean studies against the same studies for the comparator, with equal study weighting; pooled row order is descriptive. Geneformer has no clean-study estimate. All intervals condition on saved models.\n\nGSE193066's tie ceiling is retained only in conditional_GSE193066_bounds.tsv because biological identity is unresolved. No repeat-biopsy ICC or serial-biopsy label-error ceiling enters these primary panels. Other tie ceilings describe recorded label distributions, not accuracy limits or genetically verified participants.\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--apply-current-crosswalks", action="store_true")
    args = parser.parse_args()
    args.out = args.out.resolve()
    args.out.mkdir(parents=True)
    evidence = verify_evidence()
    cells, gains, counts = exposure_summary(args.out, evidence)
    framework_candidate(args.out / "framework7", cells, gains)
    shutil.copytree(TWO_MODEL, args.out / "two_model7")
    (args.out / "two_model7/REVIEW_NOTES.md").write_text("Candidate only; unchanged saved development results and all six original panels retained for separate review.\n\nThe 99 GSE267145 participant IDs are unique in the authoritative paired-assay join, but this does not establish genotype-based independence. The identity review covers neither this source nor new molecular cohorts; no automated donor exclusion is justified. Outcomes remain source-exposed development evidence. The 22 matched pairs cover fibrosis 0-1 only. No source eligibility upgrade, model fitting, weight changes, or protected Currin fold-0 rows occur in this preparation.\n\nGSE128072 remains ineligible for released-count-input validation because RNA is TPM and H3K27ac intervals use hg19. Wenz's final 175 donor crosswalk, allele mismatches and checkpoint exposure remain unresolved; no external score is claimed.\n")
    source = ROOT / "GWAS/finemapping/results/alphagenome_program/c2-endogenous-head-20260914T194000Z/tables/head_summary.tsv"
    head = pd.read_csv(source, sep="\t")
    row = head[(head.stratum == "primary_all_variants") & (head.arm == "hyenadna/delta_ridge")].iloc[0]
    rounded = "{:.4f} [{:.4f}, {:.4f}]".format(row.gain_vs_allele_identity, row.gain_ci_low, row.gain_ci_high)
    if rounded != "0.1302 [0.1161, 0.1440]":
        raise ValueError("caQTL rounding source differs")
    (args.out / "caQTL_corrected_text.md").write_text("HyenaDNA caQTL macro-Spearman gain over allele identity is **" + rounded + "** on 32,322 leads in 2,643 blocks. This corrects rounding only; the named historical C2 execution is unchanged.\n")
    changed = correct_crosswalks(args.out, evidence, args.apply_current_crosswalks)
    write_json(args.out / "preparation.json", {"status": "candidate_prepared_not_adopted", "exposure_counts": counts,
        "identity_libraries": 42, "identity_donors": 19, "models_fit": False, "weights_changed": False,
        "fold0_outcome_rows_parsed": 0, "caQTL_rounded": rounded, "modified_current_config": changed,
        "bootstrap_source_sha256": sha(EVALUATION / "bootstrap_distributions.npz"), "script_sha256": sha(Path(__file__))})
    print(json.dumps({"out": str(args.out), "counts": counts, "modified_config": len(changed)}))


if __name__ == "__main__":
    main()
