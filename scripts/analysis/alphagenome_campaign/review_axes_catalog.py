#!/usr/bin/env python3
"""Independent compute-node checks of molecular axes and actual catalog records."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sqlite3

import numpy as np
import pandas as pd

from catalog_query import query, variant_id

ROOT = Path(__file__).resolve().parents[3]
BENCH = ROOT/"Analysis/MASLD_Model_Benchmark"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--b2", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        raise SystemExit("run executable checks on a compute node")
    args.out.mkdir(parents=True, exist_ok=False)
    fix = BENCH/"executions/model-data-064-21079902/fixture"
    region_axis = pd.read_csv(fix/"molecular/h3k27ac_feature_axis.tsv", sep="\t")
    reg = pd.read_csv(args.b2/"regions.tsv", sep="\t")
    actual = region_axis.set_index("h3k27ac_feature_index").loc[reg.region_index, "opaque_source_feature_key"]
    assert list(actual) == reg.region_key.tolist(), "H3 array axis/key mismatch"
    for r in reg.itertuples():
        chrom, interval = r.region_key.split(":")
        start1, end1 = map(int, interval.split("-"))
        assert (chrom, start1-1, end1) == (r.chrom, r.start0, r.end0)
    part = pd.read_csv(fix/"molecular/participant_axis.tsv", sep="\t")
    folds = pd.read_csv(fix/"folds/participant_outer_folds.tsv", sep="\t")
    original = BENCH/"executions/corgi-film-h3k27ac-finetune-20260907T191151Z/out/film_pred_fold0.npz"
    with np.load(original, allow_pickle=False) as old:
        assert old["participant_ids"].tolist() == part.participant_id.tolist()
        assert np.array_equal(old["fold"], folds.set_index("participant_id").loc[part.participant_id, "outer_fold"].to_numpy())
        train_tiles = set(old["train_tiles"].tolist())
        held_tiles = set(old["held_tiles"].tolist())
    assert set(reg.loc[reg.region_role == "train", "tile_index"]) == train_tiles
    assert set(reg.loc[reg.region_role == "held", "tile_index"]) == held_tiles
    assert len(train_tiles) == len(held_tiles) == 32
    report = {"status":"pass", "B2_axis": {"participants":len(part), "regions":len(reg),
               "fixture_region_key_match":True,"counted_interval_conversion_match":True,
               "historical_participant_order_and_folds_match":True,"historical_32_32_tile_split_match":True}}
    dbpath = args.catalog/"regulatory.sqlite"
    with sqlite3.connect(dbpath.resolve().as_uri()+"?mode=ro", uri=True) as db:
        db.row_factory = sqlite3.Row
        counts = {name: db.execute(f"SELECT count(*) FROM {name}").fetchone()[0]
                  for name in ("source","entity","evidence","evidence_entity","relationship")}
        expected = json.loads((args.catalog/"summary.json").read_text())
        assert all(expected[k] == v for k,v in counts.items())
        assert db.execute("PRAGMA foreign_key_check").fetchall() == []
        assert db.execute("SELECT count(*) FROM entity WHERE kind='program'").fetchone()[0] == 117
        assert db.execute("SELECT count(*) FROM entity WHERE trim(identifier)='' OR identifier IS NULL").fetchone()[0] == 0, "missing molecular identities must not become shared entities"
        unsigned = db.execute("SELECT sign_convention,effect FROM evidence WHERE assay='splice_site_usage'").fetchall()
        assert unsigned and all(r["sign_convention"] == "unsigned_magnitude_no_direction" for r in unsigned)
        assert all(r["effect"] is None or r["effect"] >= 0 for r in unsigned)
        corrected = db.execute("SELECT count(*) FROM evidence WHERE assay='allelic_accessibility' AND unsupported_reason='orientation_changed_requires_reharmonization'").fetchone()[0]
        assert db.execute("SELECT count(*) FROM evidence WHERE unsupported_reason LIKE 'orientation_changed%' AND effect IS NOT NULL").fetchone()[0] == 0
        # Every available ASE record must use the same deposited allele identity
        # as its actual queried variant. No silently relocated contrast.
        ase = db.execute("SELECT e.metadata,en.identifier FROM evidence e JOIN evidence_entity ee USING(evidence_id) JOIN entity en USING(entity_key) WHERE e.assay='allelic_accessibility' AND e.effect IS NOT NULL AND en.kind='variant'").fetchall()
        for row in ase:
            inp = json.loads(row["metadata"])["input_identity"]
            assert variant_id("GRCh38:"+inp) == row["identifier"]
        zeros = db.execute("SELECT count(*) FROM evidence WHERE effect=0 AND state='recorded'").fetchone()[0]
        query_summary = {}
        for kind in ("variant","gene","region","event","program","context"):
            identifier = db.execute("SELECT identifier FROM entity WHERE kind=? ORDER BY identifier LIMIT 1", (kind,)).fetchone()[0]
            result = query(dbpath, kind, identifier, limit=2)
            assert result["entities"] and result["release_state"] == "candidate_not_adopted"
            for effect in result["evidence"]:
                assert effect["units"] and effect["source"] and effect["kind"] in ("measured", "predicted")
            query_summary[kind] = {"identity":identifier,"status":result["status"],
                                   "n_evidence":result["total_evidence"],"n_relationships":len(result["relationships"])}
            (args.out/f"example_{kind}.json").write_text(json.dumps(result, indent=2)+"\n")
        report["catalog"] = {"actual_counts":counts,"six_entry_points":query_summary,
                             "unsigned_splice_rows":len(unsigned),"corrected_ASE_masked":corrected,
                             "available_ASE_identity_checked":len(ase),"recorded_exact_zero_effects":zeros,
                             "foreign_keys_valid":True,"fixed117_preserved":True}
    report["protected_outcomes_read"] = False
    (args.out/"axis_catalog_review.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
