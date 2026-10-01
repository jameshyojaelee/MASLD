#!/usr/bin/env python3
"""Independently check imported Currin molecular identities and stored effects."""
import argparse
import json
import os
from pathlib import Path
import sqlite3

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT/"GWAS/finemapping/results/alphagenome_campaign/week1-20260915"
LABELS = ROOT/"GWAS/finemapping/results/alphagenome_program/c2-endogenous-head-20260914T194000Z/inputs/currin_lead_labels.tsv.gz"
PRED = BASE/"model/frozen/comparisons/predictions.tsv.gz"
SE = BASE/"data/endpoints-revised-21772544/caqtl_source_leads.tsv.gz"
RECIPES = [f"frozen_{length}_{pool}_{head}" for length in (2048, 16384)
           for pool in ("variant", "symmetric", "target") for head in ("ridge", "shared_mlp64")]


def connection(path):
    conn = sqlite3.connect(path.resolve().as_uri()+"?mode=ro", uri=True)
    conn.execute("PRAGMA query_only=ON")
    return conn


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--base", type=Path, default=BASE/"catalog-21772728/regulatory.sqlite")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        raise SystemExit("Run the numerical identity check on a compute node")
    args.out.mkdir(parents=True, exist_ok=False)
    labels = pd.read_csv(LABELS, sep="\t")
    labels["key"] = labels.lead_variant_id.str.replace(r"^chr", "", regex=True)
    assert len(labels) == 32322 and labels.key.is_unique
    labels = labels.set_index("key", drop=False)
    predictions = pd.read_csv(PRED, sep="\t").set_index("key")
    assert len(predictions) == 28600 and predictions.index.is_unique
    identity_columns = ["chr", "pos_hg38", "ref", "alt", "peak_id", "peak_start_hg38", "peak_stop_hg38", "heldout_fold", "block_1mb"]
    for length in (2048, 16384):
        manifest = pd.read_csv(BASE/f"model/frozen/{length}/manifest.tsv", sep="\t").set_index("key")
        assert manifest.index.is_unique and set(manifest.index) == set(labels.index)
        pd.testing.assert_frame_equal(manifest.loc[labels.index, identity_columns], labels[identity_columns], check_dtype=False)
    nominal = pd.read_csv(SE, sep="\t")
    duplicated = nominal.duplicated(["variant_id", "target_id"], keep=False)
    nominal = nominal.loc[~duplicated].set_index(["variant_id", "target_id"])
    expected_se, se_available = {}, 0
    for row in labels.itertuples(index=False):
        key = (row.lead_variant_id, row.peak_id)
        value = None
        if key in nominal.index:
            n = nominal.loc[key]
            if (np.isclose(row.beta_source, n.beta_nominal, rtol=1e-5, atol=1e-7)
                    and n.uncertainty_state == "sqrt_source_nominal_varbeta_matched_by_variant_and_peak"
                    and np.isfinite(n.varbeta) and n.varbeta > 0
                    and np.isclose(np.sqrt(n.varbeta), n.effect_se, rtol=1e-8, atol=1e-12)):
                value = float(np.sqrt(n.varbeta))
                se_available += 1
        expected_se[row.key] = value
    conn = connection(args.catalog)
    conn.execute("ATTACH DATABASE ? AS original", (args.base.resolve().as_uri()+"?mode=ro",))
    counts = {}
    for table in ("source", "entity", "evidence", "evidence_entity", "relationship"):
        assert conn.execute(f"SELECT * FROM original.{table} EXCEPT SELECT * FROM main.{table} LIMIT 1").fetchone() is None
        counts[table] = conn.execute(f"SELECT count(*) FROM main.{table}").fetchone()[0]
    assert conn.execute("SELECT count(*) FROM entity WHERE trim(identifier)='' OR identifier IS NULL").fetchone()[0] == 0
    assert conn.execute("SELECT count(*) FROM entity WHERE kind='program'").fetchone()[0] == 117
    sources = {"measured": str(LABELS.relative_to(ROOT)), "predicted": str(PRED.relative_to(ROOT))}
    imported = {"measured": 0, "predicted": 0}
    actual_se = 0
    sql = """SELECT e.effect,e.se,e.units,e.sign_convention,e.kind,e.model_version,e.metadata,
               v.identifier,r.identifier,c.identifier
             FROM evidence e
             JOIN evidence_entity ev ON e.evidence_id=ev.evidence_id AND ev.role='variant'
             JOIN entity v ON ev.entity_key=v.entity_key
             JOIN evidence_entity er ON e.evidence_id=er.evidence_id AND er.role='measured_target_region'
             JOIN entity r ON er.entity_key=r.entity_key
             JOIN evidence_entity ec ON e.evidence_id=ec.evidence_id AND ec.role='context'
             JOIN entity c ON ec.entity_key=c.entity_key
             WHERE e.source_id=?"""
    seen = set()
    for kind, source_id in sources.items():
        for effect, se, units, sign, stored_kind, model, raw_metadata, variant, region, context in conn.execute(sql, (source_id,)):
            assert stored_kind == kind and units == "source_normalized_accessibility_slope_per_ALT_dosage"
            assert context == "Currin2025:human_liver_caQTL"
            key = variant.removeprefix("GRCh38:chr")
            expected = labels.loc[key]
            assert region == "Currin2025:GRCh38:"+expected.peak_id
            metadata = json.loads(raw_metadata)
            assert metadata["heldout_fold"] == int(expected.heldout_fold)
            identity = (key, region, kind, model)
            assert identity not in seen
            seen.add(identity)
            if kind == "measured":
                assert float(effect) == float(expected.beta_alt)
                assert sign == "positive_increases_with_ALT_dosage"
                assert metadata["biological_n"] is None and metadata["source_cohort_n"] == 138
                assert metadata["effective_n_status"] == "per_variant_not_available"
                if expected_se[key] is None:
                    assert se is None
                else:
                    assert float(se) == expected_se[key]
                    actual_se += 1
            else:
                assert model in RECIPES and key in predictions.index
                assert float(effect) == float(predictions.loc[key, model])
                assert se is None and sign == "predicted_increase_with_ALT_dosage"
                assert metadata["training_folds"] == [f for f in range(5) if f != int(expected.heldout_fold)]
                assert metadata["scaling"] == "training_partition_only"
                assert metadata["evaluation"] == "fixed_recipe_single_seed_out_of_fold_development_prediction"
                assert metadata["generalization"] == "downstream_chromosome_holdout; foundation_exposure_unresolved"
                assert metadata["fitted_head_weights"] == "not_saved_not_deployed"
            imported[kind] += 1
        assert conn.execute("SELECT count(*) FROM evidence WHERE source_id=?", (source_id,)).fetchone()[0] == imported[kind]
        assert conn.execute("""SELECT count(*) FROM evidence e JOIN evidence_entity ee USING(evidence_id)
                               JOIN entity en USING(entity_key) WHERE e.source_id=? AND en.kind='gene'""", (source_id,)).fetchone()[0] == 0
    assert imported == {"measured": 32322, "predicted": 343200}
    assert actual_se == se_available
    assert conn.execute("PRAGMA foreign_key_check").fetchall() == []
    assert conn.execute("PRAGMA integrity_check").fetchone()[0] == "ok"
    conn.close()
    result = dict(status="pass", imported=imported, measurement_se_available=actual_se,
        measurement_se_unavailable=32322-actual_se, prediction_SEs=0, all_original_rows_unchanged=True,
        frozen_manifests_match_exact_variant_peak_coordinates_and_folds=True,
        all_stored_effects_match_source=True, fixed117_retained=True, no_new_gene_effect_links=True,
        source_units="normalized_accessibility_slope_per_ALT_dosage", effective_per_variant_n="unknown",
        scope="stored_development_predictions_and_marginal_source_associations; no_deployed_heads_or_external_generalization",
        protected_outcomes_read=False, counts=counts)
    (args.out/"checks.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
