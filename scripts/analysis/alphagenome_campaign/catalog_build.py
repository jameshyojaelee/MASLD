#!/usr/bin/env python3
"""Normalize an existing candidate Gene Catalog attachment into a local SQLite specimen.

No source participant matrices, protected outcomes, or canonical Catalog writes.
This is an evidence browser, not a new genetic membership or accuracy claim.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import sqlite3

from catalog_query import variant_id

ROOT = Path(__file__).resolve().parents[3]
PROG = ROOT / "GWAS/finemapping/results/alphagenome_program"
B1 = PROG / "b1-response-layer-v2-20260914T204700Z"
SPLICE = ROOT / "GWAS/finemapping/results/alphagenome_atlas/p3a-splice-direction-20260915T141415Z/tables/splice_direction_pairs.tsv"
MEMBERS = ROOT / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_membership_v2.tsv"

SCHEMA = """
PRAGMA foreign_keys=ON;
CREATE TABLE source(source_id TEXT PRIMARY KEY,path TEXT NOT NULL,sha256 TEXT NOT NULL,terms TEXT NOT NULL,exposure TEXT NOT NULL);
CREATE TABLE entity(entity_key TEXT PRIMARY KEY,kind TEXT NOT NULL,identifier TEXT NOT NULL,display_name TEXT NOT NULL,attributes TEXT NOT NULL,UNIQUE(kind,identifier));
CREATE TABLE evidence(evidence_id TEXT PRIMARY KEY,source_id TEXT NOT NULL REFERENCES source,endpoint TEXT NOT NULL,kind TEXT NOT NULL,assay TEXT NOT NULL,effect REAL,se REAL,units TEXT NOT NULL,sign_convention TEXT NOT NULL,state TEXT NOT NULL,unsupported_reason TEXT NOT NULL,model_version TEXT NOT NULL,metadata TEXT NOT NULL);
CREATE TABLE evidence_entity(evidence_id TEXT REFERENCES evidence,entity_key TEXT REFERENCES entity,role TEXT NOT NULL,PRIMARY KEY(evidence_id,entity_key,role));
CREATE TABLE relationship(from_key TEXT REFERENCES entity,to_key TEXT REFERENCES entity,link_type TEXT NOT NULL,source_id TEXT REFERENCES source,attributes TEXT NOT NULL,PRIMARY KEY(from_key,to_key,link_type,source_id));
CREATE INDEX entity_lookup ON entity(kind,display_name);
CREATE INDEX evidence_lookup ON evidence_entity(entity_key);
CREATE INDEX relationship_reverse ON relationship(to_key);
"""


def rows(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def stable(*parts):
    return hashlib.sha256(json.dumps(parts, sort_keys=True).encode()).hexdigest()[:32]


def numeric(value):
    try:
        value = float(value)
        return value if math.isfinite(value) else None
    except (ValueError, TypeError):
        return None


class Catalog:
    def __init__(self, path):
        if path.exists():
            raise FileExistsError(path)
        self.db = sqlite3.connect(path)
        self.db.executescript(SCHEMA)

    def source(self, path, terms, exposure="development_inspected; foundation_exposure_not_fully_resolved"):
        key = str(path.relative_to(ROOT))
        self.db.execute("INSERT OR IGNORE INTO source VALUES (?,?,?,?,?)", (key, key, sha(path), terms, exposure))
        return key

    def entity(self, kind, identifier, display="", **attributes):
        if not identifier or not identifier.strip():
            raise ValueError("Missing molecular identity cannot define an entity")
        key = stable(kind, identifier)
        self.db.execute("INSERT OR IGNORE INTO entity VALUES (?,?,?,?,?)", (key, kind, identifier, display or identifier, json.dumps(attributes, sort_keys=True)))
        return key

    def link(self, source, a, b, kind, **attributes):
        self.db.execute("INSERT OR IGNORE INTO relationship VALUES (?,?,?,?,?)", (a, b, kind, source, json.dumps(attributes, sort_keys=True)))

    def evidence(self, source, identity, endpoint, kind, assay, value, units, sign, objects,
                 se=None, reason="", model="", **metadata):
        effect = numeric(value)
        se = numeric(se)
        if se is not None and se < 0:
            raise ValueError("Negative standard error")
        # Different contexts/targets have different identity tuples. Repeated signal
        # membership adds a relationship to the same observation, not a replicate.
        key = stable(source, identity, endpoint)
        data = (key, source, endpoint, kind, assay, effect, se, units, sign,
                "recorded" if effect is not None else "unsupported", reason if effect is None else "",
                model, json.dumps({"release_state": "candidate_not_adopted", "uncertainty_status": "source_se" if se is not None else "not_available", **metadata}, sort_keys=True))
        old = self.db.execute("SELECT * FROM evidence WHERE evidence_id=?", (key,)).fetchone()
        if old is not None and tuple(old) != data:
            raise ValueError(f"Conflicting duplicate evidence identity: {identity} {endpoint}")
        self.db.execute("INSERT OR IGNORE INTO evidence VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)", data)
        self.db.executemany("INSERT OR IGNORE INTO evidence_entity VALUES (?,?,?)", [(key, e, role) for e, role in objects if e])


def load_layer(catalog):
    source = catalog.source(B1 / "response_layer_variants.tsv.gz", "local_candidate; aggregate_source_conditions_apply; hosted_outputs_not_training_features")
    refresh_path = PROG / "b1-dbsnp-refresh-20260915T084700Z/tables/research_layer_orientation_refreshed_all_rows.tsv"
    refresh_source = catalog.source(refresh_path, "local_candidate_orientation_reference")
    orientations = {(r["signal_uid"], r["variant_key"]): r for r in rows(refresh_path)}
    panel = {r["exported_column"]: r for r in rows(B1 / "track_panel_composition.tsv")}
    contexts = {}
    for column, record in panel.items():
        contexts[column] = catalog.entity("context", "Atlas:" + column, column, **record,
                                          observed_or_synthetic="pretrained_model_track_panel",
                                          capability="static_allele_prediction; donor_or_treatment_interaction_unvalidated")
    for row in rows(B1 / "response_layer_variants.tsv.gz"):
        orientation = orientations.get((row["signal_uid"], row["variant_key"]), {})
        original = row["variant_uid"]
        corrected = orientation.get("resolved_identity__current_20260914T210732Z") or original
        valid = True
        try:
            identity = variant_id("GRCh38:" + corrected)
        except ValueError:
            valid = False
            identity = "unresolved:" + (row["variant_key"] or original or row["source_variant_id"])
        variant = catalog.entity("variant", identity, original_source_id=row["source_variant_id"],
            deposited_identity=original, orientation_source=orientation.get("orientation_source__current_20260914T210732Z", "") or row["mapping_status"],
            reference_verification="source_mapped_not_rechecked_in_this_build" if valid else "unresolved",
            normalized="not_independently_left_normalized", build="GRCh38", coordinate_convention="1_based_variant_position")
        signal = catalog.entity("signal", row["signal_uid"], trait=row["trait"], trait_class=row["trait_class"],
                                posterior_definition=row["posterior_definition"], covered_share=row["signal_queried_share"])
        gene = catalog.entity("gene", row["ensembl"], row["gene"], annotation_version="source_version_unresolved") if row["ensembl"] else None
        catalog.link(source, variant, signal, "signal_membership", posterior=row["weight"], posterior_definition=row["posterior_definition"])
        if corrected != original:
            catalog.link(refresh_source, variant, signal, "orientation_amendment", deposited=original, corrected=corrected)
        if gene:
            catalog.link(source, variant, gene, "candidate_signal_target", signal=row["signal_uid"], claim="colocalized_candidate_not_causation")
        objects = [(variant, "variant"), (signal, "signal"), (gene, "candidate_target")]
        regions = []
        for peak in filter(None, row["measured_peak"].split(";")):
            region = catalog.entity("region", "GRCh38:" + peak, coordinate_convention="source_peak_key; inspect_source_before_sequence_extraction")
            catalog.link(source, variant, region, "measured_peak_overlap", lineages=row["measured_lineages"])
            regions.append((region, "measured_region"))
        if row["rna_strongest_gene"]:
            alternative = catalog.entity("gene", row["rna_strongest_gene"], annotation_version="source_version_unresolved")
            catalog.link(source, variant, alternative, "strongest_predicted_RNA_target", claim="prediction_not_target_validation")
        for column, assay, sign in [
            ("atac_primary_liver_quantile", "ATAC", "signed_ALT_relative_to_REF_quantile"),
            ("dnase_primary_liver_quantile", "DNase", "signed_ALT_relative_to_REF_quantile"),
            ("h3k27ac_primary_liver_quantile", "H3K27ac", "signed_ALT_relative_to_REF_quantile"),
            ("rna_target_gene_primary_liver_quantile", "RNA", "signed_ALT_relative_to_REF_quantile"),
            ("splice_site_usage_primary_liver_quantile", "splice_site_usage", "unsigned_magnitude_no_direction")]:
            context = contexts.get(column) or catalog.entity("context", "Atlas:"+column, life_stage="unresolved")
            changed = corrected != original
            catalog.evidence(source, (row["variant_key"], row["ensembl"] if assay in {"RNA", "splice_site_usage"} else "", column),
                column, "predicted", assay, None if changed else row[column], "Atlas_quantile_dimensionless", sign,
                objects + [(context, "context")], reason="orientation_changed_requires_rescoring" if changed else "unserved_or_unavailable_track",
                model="AlphaGenome_Atlas_archived_2026-09", input_identity=original,
                applicability="not_a_donor_or_treatment_prediction")
        for accession in ("gse281367", "gse244832"):
            context = catalog.entity("context", accession + ":liver_allelic_accessibility", tissue="liver", life_stage="source_registry", observed_or_synthetic="measured")
            column = f"ase_{accession}_mean_log2_alt_over_ref"
            catalog.evidence(source, (row["variant_key"], accession), column, "measured", "allelic_accessibility",
                None if corrected != original else row[column], "log2_ALT_over_REF_mean_over_heterozygous_donors", row["ase_sign_refers_to"],
                objects+regions+[(context,"context")], se=row[f"ase_{accession}_se_log2"],
                reason="orientation_changed_requires_reharmonization" if corrected != original else "no_eligible_allelic_measurement", biological_unit="donor",
                deposited_effect_before_orientation_review=numeric(row[column]) if corrected != original else None,
                biological_n=numeric(row[f"ase_{accession}_n_het_donors"]), input_identity=original)
        # QTL signs retain the source effect allele, even when hg38 REF/ALT changed.
        if gene:
            context = catalog.entity("context", "source_liver_eQTL", observed_or_synthetic="measured", life_stage="source_specific_unresolved")
            catalog.evidence(source, (row["variant_key"], row["ensembl"], row["gwas_name"]),
                "liver_eqtl", "measured", "eQTL", row["eqtl_beta_allele1"], "source_eQTL_beta",
                row["direction_sign_refers_to"], objects+[(context,"context")], reason="no_harmonized_eQTL_effect",
                biological_unit="source_participants", biological_n=None,
                estimand="marginal_association_not_causal_effect", effect_allele_identity=row["source_variant_id"])


def load_events(catalog):
    source = catalog.source(SPLICE, "local_candidate; hosted_outputs_not_training_features")
    context = catalog.entity("context", "signed_splice:historical_liver_panel", life_stage="mixed_or_unresolved", capability="weak_direction_below_original_ranking_criteria")
    for row in rows(SPLICE):
        variant = catalog.entity("variant", variant_id("GRCh38:" + row["variant_uid"]), build="GRCh38", reference_verification="source_mapped_not_rechecked_in_this_build")
        gene = catalog.entity("gene", row["gene"], annotation_version="phenotype_ID_retains_version") if row["gene"].strip() else None
        event = catalog.entity("event", row["phenotype_id"], strand="not_present_in_archived_summary",
            intron_start=row["intron_start"], intron_end=row["intron_end"], source_convention="GTEx_LeafCutter_closed_intron",
            event_definition="historical_loose_shared_donor_or_acceptor_denominator",
            gained_lost_junctions="not_reconstructible_from_summary_new_inference_required") if row["phenotype_id"].strip() else None
        if event and gene:
            catalog.link(source, event, gene, "source_event_gene")
        objects = [(variant,"variant"),(gene,"event_gene"),(event,"event"),(context,"context")]
        catalog.evidence(source, (row["variant_uid"],row["phenotype_id"]), "excision_ratio", "predicted", "splicing",
            row["excision_ratio_log2"], "log2_ALT_over_REF_excision_ratio", "ALT_relative_to_REF",
            objects, reason=row["state"] if event else "event_identity_missing", model="AlphaGenome_model_API_archived_2026-09-15",
            source_state=row["state"], event_identity_status="identified" if event else "missing",
            historical_decision="weak_direction_below_ranking_criteria; magnitude_and_event_detection_failed",
            n_cluster=numeric(row["n_cluster"]))
        catalog.evidence(source, (row["variant_uid"],row["phenotype_id"]), "source_sqtl", "measured", "sQTL",
            row["slope"], "source_sQTL_slope", "source_effect_allele_orientation_see_source",
            objects, reason="missing_sQTL_effect" if event else "event_identity_missing", biological_unit="source_participants", biological_n=None)


def load_programs(catalog):
    source = catalog.source(MEMBERS, "fixed117_membership_unchanged; no_outcome_optimization", "fixed_before_campaign")
    ids = set()
    for row in rows(MEMBERS):
        ids.add(row["program_uid"])
        program = catalog.entity("program", row["program_uid"], discovery_lineage=row["cell_type"], module=row["module"], membership_sha256=row["membership_sha256"])
        gene = catalog.entity("gene", "symbol:" + row["canonical_gene"], row["canonical_gene"], annotation_version="fixed117_source_symbol; no_unverified_Ensembl_join")
        catalog.link(source, program, gene, "fixed_program_member", original_l1_weight=row["original_l1_weight"], canonical_weight_text=row["canonical_weight_text"])
    if len(ids) != 117:
        raise ValueError(f"Fixed program count changed: {len(ids)}")


def build(out):
    out.mkdir(parents=True, exist_ok=True)
    catalog = Catalog(out / "regulatory.sqlite")
    try:
        load_layer(catalog)
        load_events(catalog)
        load_programs(catalog)
        catalog.db.commit()
        if catalog.db.execute("PRAGMA foreign_key_check").fetchall():
            raise ValueError("Broken scientific identity relationship")
        counts = {}
        for table in ("source", "entity", "evidence", "evidence_entity", "relationship"):
            cursor = catalog.db.execute(f"SELECT * FROM {table}")
            with gzip.open(out / (table + ".tsv.gz"), "wt") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow([c[0] for c in cursor.description])
                writer.writerows(cursor)
            counts[table] = catalog.db.execute(f"SELECT count(*) FROM {table}").fetchone()[0]
        counts["by_entity"] = dict(catalog.db.execute("SELECT kind,count(*) FROM entity GROUP BY kind"))
        counts.update(release_state="candidate_not_adopted", protected_outcomes_read=False,
                      source_data_rehosted=False, joint_haplotype_inference="separate_local_SNV_route; not_reconstructed_from_catalog_records",
                      scope="bounded_existing_aggregate_attachment_not_full_Atlas")
        (out / "summary.json").write_text(json.dumps(counts, indent=2) + "\n")
        (out / "schema.sql").write_text(SCHEMA)
        print(json.dumps(counts, indent=2))
    finally:
        catalog.db.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    build(parser.parse_args().out)
