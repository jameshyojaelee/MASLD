#!/usr/bin/env python3
"""Add checked Currin measurements and development predictions to a new specimen.

All twelve fixed recipes are retained. This imports stored out-of-fold values;
it does not select a model, deploy fitted heads, or create genetic gene links.
"""
import argparse
import csv
import gzip
import json
import os
from pathlib import Path
import platform
import sqlite3

import numpy as np
import pandas as pd

from catalog_build import Catalog, ROOT, SCHEMA, numeric, sha
from catalog_query import connect, query, variant_id

BASE = ROOT / "GWAS/finemapping/results/alphagenome_campaign/week1-20260915"
LABELS = ROOT / "GWAS/finemapping/results/alphagenome_program/c2-endogenous-head-20260914T194000Z/inputs/currin_lead_labels.tsv.gz"
SOURCE_SE = BASE / "data/endpoints-revised-21772544/caqtl_source_leads.tsv.gz"
PREDICTIONS = BASE / "model/frozen/comparisons/predictions.tsv.gz"
CONFIGS = BASE / "model/configs"
RECIPES = [f"frozen_{length}_{pool}_{head}" for length in (2048, 16384)
           for pool in ("variant", "symmetric", "target") for head in ("ridge", "shared_mlp64")]
TABLES = ("source", "entity", "evidence", "evidence_entity", "relationship")
UNITS = "source_normalized_accessibility_slope_per_ALT_dosage"


def matched_uncertainty(labels, source):
    """Retain unavailable uncertainty without substituting another target."""
    columns = ["variant_id", "target_id", "beta_nominal", "varbeta", "effect_se", "uncertainty_state"]
    source = source[columns].copy()
    duplicates = source.duplicated(["variant_id", "target_id"], keep=False)
    unique = source.loc[~duplicates].rename(columns={"variant_id": "lead_variant_id", "target_id": "peak_id"})
    joined = labels.merge(unique, on=["lead_variant_id", "peak_id"], how="left", validate="one_to_one", indicator=True)
    positive = np.isfinite(joined.varbeta) & joined.varbeta.gt(0)
    rooted = np.sqrt(joined.varbeta.where(positive))
    status = np.select([
        joined._merge.ne("both"),
        ~np.isclose(joined.beta_source, joined.beta_nominal, rtol=1e-5, atol=1e-7),
        joined.uncertainty_state.ne("sqrt_source_nominal_varbeta_matched_by_variant_and_peak"),
        ~positive,
        ~np.isclose(rooted, joined.effect_se, rtol=1e-8, atol=1e-12)],
        ["missing_or_nonunique_exact_variant_peak", "source_beta_mismatch",
         "uncertainty_definition_unverified", "nonpositive_or_missing_varbeta",
         "deposited_SE_disagrees_with_varbeta"], default="verified_exact_variant_peak_native_varbeta")
    joined["verified_se"] = rooted.where(status == "verified_exact_variant_peak_native_varbeta")
    joined["se_status"] = status
    return joined.drop(columns="_merge")


def validate_predictions(labels, predictions):
    if not labels.lead_variant_id.is_unique or not predictions.key.is_unique:
        raise ValueError("Duplicate variant identities cannot be treated as independent predictions")
    expected = {"key", "block_1mb", "heldout_fold", "beta_alt", *RECIPES}
    if set(predictions.columns) != expected:
        raise ValueError("The complete twelve-recipe comparison is required")
    source = labels.assign(key=labels.lead_variant_id.str.removeprefix("chr")).set_index("key")
    source = source.loc[predictions.key]
    np.testing.assert_allclose(predictions.beta_alt, source.beta_alt, rtol=1e-6, atol=1e-7)
    np.testing.assert_array_equal(predictions.heldout_fold, source.heldout_fold)
    np.testing.assert_array_equal(predictions.block_1mb, source.block_1mb)
    if not np.isfinite(predictions[["beta_alt", *RECIPES]].to_numpy(dtype=float)).all():
        raise ValueError("Incomplete predictions must not become recorded zero effects")
    if set(predictions.heldout_fold) != set(range(5)):
        raise ValueError("Expected five established development folds")


def copy_specimen(base, destination):
    if destination.exists():
        raise FileExistsError(destination)
    catalog = Catalog.__new__(Catalog)
    catalog.db = sqlite3.connect(destination)
    with connect(base) as original:
        original.backup(catalog.db)
    catalog.db.execute("PRAGMA foreign_keys=ON")
    return catalog


def write_exports(catalog, out):
    counts = {}
    for table in TABLES:
        cursor = catalog.db.execute(f"SELECT * FROM {table}")
        with gzip.open(out / (table + ".tsv.gz"), "wt") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow([c[0] for c in cursor.description])
            writer.writerows(cursor)
        counts[table] = catalog.db.execute(f"SELECT count(*) FROM {table}").fetchone()[0]
    counts["by_entity"] = dict(catalog.db.execute("SELECT kind,count(*) FROM entity GROUP BY kind"))
    return counts


def build(args):
    if not os.environ.get("SLURM_JOB_ID"):
        raise SystemExit("Run the build and scientific checks on a compute node")
    args.out.mkdir(parents=True, exist_ok=False)
    base_hash = sha(args.base)
    labels = pd.read_csv(LABELS, sep="\t")
    predictions = pd.read_csv(PREDICTIONS, sep="\t")
    if len(labels) != 32322 or len(predictions) != 28600:
        raise ValueError("Unexpected eligible or matched denominator")
    validate_predictions(labels, predictions)
    np.testing.assert_allclose(labels.beta_alt, labels.beta_source, rtol=0, atol=0)
    if not labels.reference_match.eq(True).all() or not labels.window_ok.eq(True).all():
        raise ValueError("Source reference or sequence eligibility has changed")
    source = pd.read_csv(SOURCE_SE, sep="\t")
    labels = matched_uncertainty(labels, source)
    config_paths = [CONFIGS / (name + ".json") for name in RECIPES]
    registry = {}
    for name, path in zip(RECIPES, config_paths):
        recipe = json.loads(path.read_text())
        if recipe["id"] != name or recipe["mode"] != "frozen":
            raise ValueError("Wrong recipe definition")
        meta_path = BASE / f"model/frozen/{recipe['length']}/complete.json"
        meta = json.loads(meta_path.read_text())
        if not meta["complete"] or meta["allele_order"] != ["REF", "ALT"]:
            raise ValueError("Unverified sequence representation")
        registry[name] = dict(recipe=recipe, recipe_sha256=sha(path), representation=meta,
            representation_receipt_sha256=sha(meta_path), fitted_head_weights="not_saved_not_deployed",
            prediction_source_sha256=sha(PREDICTIONS),
            fit_code_sha256=sha(Path(__file__).with_name("model_frozen_fit.py")))
    catalog = copy_specimen(args.base, args.out / "regulatory.sqlite")
    try:
        before = {table: catalog.db.execute(f"SELECT count(*) FROM {table}").fetchone()[0] for table in TABLES}
        measured_source = catalog.source(LABELS, "public_aggregate_association_source_conditions_apply")
        uncertainty_source = catalog.source(SOURCE_SE, "source_nominal_varbeta_exact_variant_peak_join")
        predicted_source = catalog.source(PREDICTIONS, "AlphaGenome_derived_noncommercial_terms; aggregate_development_predictions")
        context = catalog.entity("context", "Currin2025:human_liver_caQTL", "Currin 2025 liver accessibility",
            tissue="human_liver", assay="ATAC_caQTL", biological_unit="source_participant",
            source_cohort_participants=138, effective_per_variant_n="not_available",
            applicability="source_liver_association; not_patient_MASLD_response",
            population="historical_q_lt_0.05_autosomal_SNV_leads_selected_by_source; one_target_per_variant",
            competing_explanations="marginal_association_can_reflect_linked_variants; no_causal_target_gene_established")
        objects = {}
        for row in labels.itertuples(index=False):
            identifier = variant_id("GRCh38:" + row.lead_variant_id)
            variant = catalog.entity("variant", identifier, build="GRCh38",
                coordinate_convention="1_based_variant_position", reference_verification="Currin_C2_reference_checked_SNV",
                normalized="biallelic_SNV")
            region = catalog.entity("region", "Currin2025:GRCh38:" + row.peak_id, row.peak_id,
                build="GRCh38", chromosome=row.chr, source_start=int(row.peak_start_hg38),
                source_end=int(row.peak_stop_hg38),
                coordinate_convention="source_BED_peak_coordinates_retained; no_new_sequence_extraction",
                molecular_role="measured_accessibility_target_not_gene_link")
            obj = [(variant, "variant"), (region, "measured_target_region"), (context, "context")]
            key = row.lead_variant_id.removeprefix("chr")
            objects[key] = (obj, row.peak_id)
            catalog.link(measured_source, variant, region, "source_marginal_caQTL_association",
                peak_id=row.peak_id, claim="association_not_causal_regulation")
            catalog.evidence(measured_source, (identifier, row.peak_id), "accessibility_caQTL", "measured", "caQTL",
                row.beta_alt, UNITS, "positive_increases_with_ALT_dosage", obj, se=row.verified_se,
                estimand="marginal_source_association", biological_unit="source_participant", biological_n=None,
                source_cohort_n=138, effective_n_status="per_variant_not_available",
                source_q=numeric(row.q_val), nominal_p=numeric(row.p_nominal),
                maf=numeric(row.maf), imputation_r2=numeric(row.imputation_r2),
                heldout_fold=int(row.heldout_fold), historical_1mb_bin=row.block_1mb,
                uncertainty_status=row.se_status, uncertainty_source=uncertainty_source,
                population="existing_significance_selected_leads; not_all_nominal_variant_peak_pairs")
        for name in RECIPES:
            recipe = registry[name]["recipe"]
            for row in predictions.itertuples(index=False):
                obj, peak = objects[row.key]
                catalog.evidence(predicted_source, (row.key, peak, name), "accessibility_caQTL", "predicted", "caQTL",
                    getattr(row, name), UNITS, "predicted_increase_with_ALT_dosage", obj, model=name,
                    evaluation="fixed_recipe_single_seed_out_of_fold_development_prediction",
                    heldout_fold=int(row.heldout_fold), training_folds=[f for f in range(5) if f != row.heldout_fold],
                    scaling="training_partition_only", sequence_length_bp=recipe["length"], pooling=recipe["pooling"],
                    head=recipe["head"], seed=recipe["seed"], recipe_sha256=registry[name]["recipe_sha256"],
                    uncertainty_status="individual_prediction_uncertainty_not_estimated; measurement_SE_not_reused",
                    model_registry="model_registry.json", fitted_head_weights="not_saved_not_deployed",
                    generalization="downstream_chromosome_holdout; foundation_exposure_unresolved",
                    population="28600_identical_rows_across_all_twelve_recipes; coverage_not_full_population")
            print(json.dumps(dict(imported_recipe=name, rows=len(predictions))), flush=True)
        catalog.db.commit()
        if catalog.db.execute("PRAGMA foreign_key_check").fetchall():
            raise ValueError("Broken molecular identity relationships")
        if catalog.db.execute("PRAGMA integrity_check").fetchone()[0] != "ok":
            raise ValueError("SQLite integrity check failed")
        counts = write_exports(catalog, args.out)
        assert counts["evidence"] - before["evidence"] == 32322 + 12 * 28600
        assert counts["by_entity"]["program"] == 117
        assert sha(args.base) == base_hash, "Original specimen changed during import"
        registry_path = args.out / "model_registry.json"
        registry_path.write_text(json.dumps(registry, indent=2) + "\n")
        (args.out / "schema.sql").write_text(SCHEMA)
        summary = dict(**counts, previous_counts=before, imported_measurements=32322,
            imported_predictions=12*28600, matched_variants=28600, recipes=RECIPES,
            measurement_SE_status=labels.se_status.value_counts().to_dict(),
            release_state="candidate_not_adopted", protected_outcomes_read=False,
            source_participant_data_rehosted=False, head_weights_deployed=False,
            original_specimen=str(args.base), original_sha256=base_hash,
            source_hashes={str(p): sha(p) for p in (LABELS, SOURCE_SE, PREDICTIONS, registry_path)},
            environment=dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__, sqlite=sqlite3.sqlite_version))
        (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    finally:
        catalog.db.close()
    # Outcome-independent first deposited matched row, not the largest effect.
    example_key = "GRCh38:chr" + predictions.iloc[0].key
    result = query(args.out / "regulatory.sqlite", "variant", example_key, assay="caQTL")
    (args.out / "example_currin_variant.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, default=BASE / "catalog-21772728/regulatory.sqlite")
    parser.add_argument("--out", type=Path, required=True)
    build(parser.parse_args())
