#!/usr/bin/env python3
"""Correct local RNA eligibility and expose both COLOC methods without reclassifying genes."""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import html
import json
from pathlib import Path
import sys

import pandas as pd


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, records):
    with path.open("x") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root, out = args.project_root.resolve(), args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    sys.path.insert(0, str(root / "scripts/portal"))
    import build_evidence_passport_bundle as adapter
    import build_resource_catalog_sidecar as sidecar

    source_dir = root / "figures/candidates/catalog-r8-completion-20261014T141645Z_22236527"
    original = source_dir / "molecular_objects/molecular_object_sidecar.jsonl"
    bundle = source_dir / "catalog_local"
    index_path = bundle / "passport_gene_index.parquet"
    evidence_path = bundle / "passport_evidence_long.parquet"
    adjud_path = root / "RNA-seq/results/evidence_passports/candidates/gene-catalog-corrected-r8-20261009T142027Z/preflight/gen_identity_adjudication_v1/gen_frozen_classes_ensembl_adjudicated.tsv"
    native_path = root / "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv"
    registry_path = root / "GWAS/finemapping/config/gwas_registry.tsv"
    genes = pd.read_parquet(index_path)
    evidence = pd.read_parquet(evidence_path)
    assert len(genes) == 28233 and len(evidence) == 56466
    assert genes.ensembl_id.is_unique
    adjud = {r["adjudicated_ensembl_id"]: r for r in read_tsv(adjud_path)}
    registry = {r["study_name"]: r for r in read_tsv(registry_path)}
    with native_path.open() as handle:
        native_rows = list(csv.DictReader(handle))
    native = {r["ensembl"].split(".")[0]: r for r in native_rows}
    assert len(native) == len(native_rows)
    ev = {(str(r.ensembl_id), str(r.evidence_domain)): r._asdict() for r in evidence.itertuples(index=False)}
    assert len(ev) == len(evidence)

    changed = []
    old_classes, new_classes = Counter(), Counter()
    old_ids, new_ids, methods, methods_by_id, identity_rows, examples = [], [], [], {}, [], []
    with original.open() as src, (out / "molecular_object_sidecar.jsonl").open("x") as dst:
        for line in src:
            before = json.loads(line)
            after = json.loads(line)
            ens = before["ensembl_id_base"]
            rna = ev[ens, "transcriptomics"]
            eligibility = json.loads(after["assay_applicability"])
            expected = "untestable" if rna["call_state"] == "untestable" or rna["testability_state"] in {"untestable", "not_detected"} else "indeterminate"
            # Preserve annotation-held or reference-excluded dispositions.
            if after["counting_annotation_disposition"] == "held_untestable":
                expected = "untestable"
            elif after["counting_annotation_disposition"] == "excluded_reference":
                expected = "not_applicable"
            if eligibility["bulk_rna"] != expected:
                changed.append({"ensembl_id": ens, "symbol": before["symbol"], "rna_call": rna["call_state"],
                                "input_testability_state": rna["testability_state"], "old_applicability": eligibility["bulk_rna"],
                                "corrected_applicability": expected})
            eligibility["bulk_rna"] = expected
            after["assay_applicability"] = sidecar.catalog.canonical_assay_applicability(eligibility)
            # No gene call, membership, route, record identity or native measurement changes.
            untouched = set(before) - {"assay_applicability"}
            assert all(before[k] == after[k] for k in untouched)
            old_classes[before["primary_evidence_class"]] += 1
            new_classes[after["primary_evidence_class"]] += 1
            old_ids.append(before["catalog_object_id"])
            new_ids.append(after["catalog_object_id"])
            dst.write(json.dumps(after, sort_keys=True, allow_nan=False) + "\n")

            source = adjud[ens]
            raw = native.get(ens, {})
            susie = adapter._optional_float(raw.get("coloc_best_susie_pp4"), "susie")
            abf = adapter._optional_float(raw.get("coloc_best_abf_pp4"), "abf")
            susie_study = raw.get("susie_driving_gwas", "") if susie is not None else ""
            abf_study = raw.get("abf_driving_gwas", "") if abf is not None else ""
            assert not (susie is not None and not susie_study)
            assert not (abf is not None and not abf_study)
            for study in (susie_study, abf_study):
                assert not study or study in registry, study
            display_susie = adapter._optional_float(source["coloc_best_susie_pp4"], "display_susie")
            display_abf = adapter._optional_float(source["coloc_best_abf_pp4"], "display_abf")
            genetics = ev[ens, "genetics"]
            assay_method = "SuSiE_COLOC" if genetics["call_state"] == "untestable" else adapter._genetic_display_method(source, display_susie, display_abf)
            methods_by_id[ens] = assay_method
            display_method = assay_method if pd.notna(genetics["estimate"]) else "not_available"
            row = {"ensembl_id": ens, "symbol": before["symbol"], "primary_evidence_class": before["primary_evidence_class"],
                   "r8_genetics_call_state": genetics["call_state"], "class_call_method": "Tier-1/2 SuSiE COLOC only",
                   "native_stable_id_result": "available" if raw else "no_same_stable_id_result",
                   "source_adjudication_rule": source["identity_resolution_rule"],
                   "r8_display_method": display_method, "r8_display_posterior": genetics["estimate"] if pd.notna(genetics["estimate"]) else None,
                   "r8_display_study": genetics["dataset_id"] if pd.notna(genetics["estimate"]) else "",
                   "susie_pp_h4": susie, "susie_study": susie_study,
                   "susie_ancestry": registry[susie_study]["ancestry"] if susie_study else "",
                   "abf_pp_h4": abf, "abf_study": abf_study,
                   "abf_ancestry": registry[abf_study]["ancestry"] if abf_study else "",
                   "gen_susie_state_t12": raw.get("susie_state_t12", "no_same_stable_id_result"),
                   "limitation": "Method-specific results join by native stable Ensembl ID. ABF is descriptive and does not create an r8 genetic support call. No class, route or membership changes."}
            methods.append(row)
            identity_rows.append({"ensembl_id": ens, "symbol": before["symbol"], "primary_evidence_class": before["primary_evidence_class"],
                                  "r8_genetics_call_state": genetics["call_state"], "gen_susie_state_t12": row["gen_susie_state_t12"],
                                  "gen_coloc_best_susie_pp4": susie})
            if before["symbol"] in {"GNMT", "GLP1R", "THRB", "HKDC1", "MTARC1"}:
                examples.append({"record": after, "methods": row})
    assert len(old_ids) == len(set(old_ids)) == 28233 and old_ids == new_ids
    assert old_classes == new_classes
    assert len(changed) == 4218
    assert all(r["rna_call"] == "untestable" and r["corrected_applicability"] == "untestable" for r in changed)
    glp1r = next(x for x in examples if x["record"]["symbol"] == "GLP1R")
    assert json.loads(glp1r["record"]["assay_applicability"])["bulk_rna"] == "untestable"
    gnmt = next(x for x in examples if x["record"]["symbol"] == "GNMT")["methods"]
    assert gnmt["r8_display_method"] == "ABF_COLOC" and gnmt["susie_pp_h4"] < .5 < gnmt["abf_pp_h4"]
    corrected_evidence = evidence.copy()
    genetic_mask = corrected_evidence.evidence_domain == "genetics"
    corrected_evidence.loc[genetic_mask, "assay"] = corrected_evidence.loc[genetic_mask, "ensembl_id"].map(methods_by_id)
    pd.testing.assert_frame_equal(evidence.drop(columns="assay"), corrected_evidence.drop(columns="assay"))
    corrected_evidence.to_parquet(out / "passport_evidence_long_method_labeled.parquet", index=False)
    write_tsv(out / "rna_eligibility_changes.tsv", changed)
    write_tsv(out / "genetics_method_evidence.tsv", methods)
    write_tsv(out / "catalog_identity_current_classes.tsv", identity_rows)
    # Verify the changed producer itself, including the not-detected vocabulary mismatch.
    probe = {"primary_evidence_class": "untested", "ensembl_id": "ENSG00000112164", "symbol": "GLP1R", "next_experiment_rule_id": "EXP_DIRECT_ASSAY_V1", "limitation": ""}
    probe_out = sidecar.object_record(probe, None, None, {"transcriptomics": {"call_state": "untestable", "testability_state": "not_detected"}}, "FIXTURE", True)
    assert json.loads(probe_out["assay_applicability"])["bulk_rna"] == "untestable"
    for row in methods:
        if row["r8_genetics_call_state"] == "supported":
            assert row["r8_display_method"] == "SuSiE_COLOC"
    checks = {"status": "PASS", "gene_records": len(old_ids), "gene_id_order_unchanged": old_ids == new_ids,
              "classes_unchanged": old_classes == new_classes, "class_counts": dict(new_classes),
              "rna_applicability_corrected": len(changed), "all_other_companion_fields_unchanged": True,
              "all_evidence_values_classes_and_routes_unchanged": True,
              "assay_context_method_counts": dict(Counter(methods_by_id.values())),
              "displayed_posterior_method_counts": dict(Counter(r["r8_display_method"] for r in methods)),
              "native_method_values_do_not_use_symbol_only_joins": True,
              "GLP1R_bulk_rna": "untestable", "GNMT_methods": gnmt,
              "source_inputs": {str(p): digest(p) for p in (original, index_path, evidence_path, adjud_path, native_path, registry_path)},
              "producer_files": {str(root / p): digest(root / p) for p in ("scripts/portal/build_resource_catalog_sidecar.py", "scripts/portal/build_evidence_passport_bundle.py")}}
    (out / "checks.json").write_text(json.dumps(checks, indent=2, allow_nan=False) + "\n")
    (out / "worked_entries.json").write_text(json.dumps(examples, indent=2, allow_nan=False) + "\n")
    rows = "".join("<tr>" + "".join("<td>" + html.escape(str(x["methods"][k])) + "</td>" for k in ("symbol", "primary_evidence_class", "susie_pp_h4", "susie_study", "abf_pp_h4", "abf_study")) + "</tr>" for x in examples)
    (out / "index.html").write_text('<meta charset="utf-8"><title>MASLD Catalog corrected companion</title><style>body{font:16px sans-serif;margin:2rem;max-width:1100px}td,th{padding:.5rem;text-align:left;border-bottom:1px solid #ddd}</style><h1>MASLD Catalog corrected companion</h1><p>28,233 identities and all r8 classes, calls, routes and memberships preserved. Bulk RNA applicability corrected for 4,218 untestable genes. SuSiE and ABF posteriors are separate; only SuSiE determines the r8 genetics support call. Posterior estimates refer to their named studies and ancestry in the downloadable table.</p><table><tr><th>Gene</th><th>r8 class</th><th>SuSiE PP.H4</th><th>SuSiE study</th><th>ABF PP.H4</th><th>ABF study</th></tr>' + rows + '</table><ul><li><a href="molecular_object_sidecar.jsonl">Corrected complete molecular-object records</a></li><li><a href="genetics_method_evidence.tsv">Method-specific genetics for every Catalog identity</a></li><li><a href="passport_evidence_long_method_labeled.parquet">Separate assay evidence with explicit method labels</a></li><li><a href="rna_eligibility_changes.tsv">RNA applicability corrections</a></li><li><a href="checks.json">Checks and source files</a></li></ul><p>Local review candidate. Existing named outputs are preserved. No gene class or causal mechanism is inferred from ABF reference values.</p>')
    print(json.dumps({k: checks[k] for k in ("status", "gene_records", "class_counts", "rna_applicability_corrected", "displayed_posterior_method_counts")}))


if __name__ == "__main__":
    main()
