#!/usr/bin/env python3
"""Build the adopted priority-eligibility Catalog from preserved source tables."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import html
import json
from pathlib import Path
import shutil

import pandas as pd

import build_evidence_passport_ui as ui
import generate_evidence_passports as contract
import gene_catalog_v2 as molecular


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def tsv(path):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def save_tsv(frame, path):
    frame.to_csv(path, sep="\t", index=False, lineterminator="\n")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[2])
    p.add_argument("--config", type=Path, default=Path(__file__).with_name("catalog_current.json"))
    p.add_argument("--output", type=Path)
    a = p.parse_args()
    root = a.project_root.resolve()
    c = json.loads(a.config.read_text())
    assert c["adoption_status"] == "adopted" and c["expected_class_changes"] == 180
    out = (a.output or root / c["output"]).resolve()
    out.mkdir(parents=True, exist_ok=False)
    rel = c["release_id"]
    source = root / c["source_bundle"]
    companion = root / c["corrected_companion"]
    annotation = root / c["annotation_bundle"]
    inputs = {"config": a.config.resolve(), "proposal": root / c["eligibility_proposals"], "source_index": source / "passport_gene_index.parquet",
              "source_evidence": source / "passport_evidence_long.parquet", "source_coverage": source / "passport_coverage_long.parquet",
              "source_membership": source / "passport_program_membership.parquet", "corrected_molecular": companion / "molecular_object_sidecar.jsonl",
              "annotation_r9d": annotation / "catalog_r9d_gene_table.parquet"}
    before_hash = {k: sha(path) for k, path in inputs.items()}
    genes = pd.read_parquet(inputs["source_index"])
    old_genes = genes.copy(deep=True)
    evidence = pd.read_parquet(companion / "passport_evidence_long_method_labeled.parquet")
    old_evidence = pd.read_parquet(inputs["source_evidence"])
    coverage = pd.read_parquet(inputs["source_coverage"])
    proposals = tsv(inputs["proposal"])
    moves = proposals[proposals.priority_rule_eligibility_transition == "True"].set_index("ensembl")
    assert len(moves) == c["expected_genetics_state_changes"] and moves.index.is_unique
    old_genetics = old_evidence[old_evidence.evidence_domain == "genetics"].set_index("ensembl_id")
    for ens, move in moves.iterrows():
        gene = old_genes[old_genes.ensembl_id == ens].iloc[0]
        assert gene.primary_evidence_class == move.named_catalog_class
        assert old_genetics.loc[ens, "call_state"] == move.named_catalog_genetic_evidence_state == "indeterminate"
    changed_gene_ids = set(moves.index[moves.priority_rule_class_would_change == "True"])
    assert len(changed_gene_ids) == c["expected_class_changes"]
    gene_mask = genes.ensembl_id.isin(changed_gene_ids)
    genes.loc[gene_mask, "primary_evidence_class"] = "untested"
    for i in genes.index[gene_mask]:
        role, role_rule, experiment = contract.route_hypothesis("untested")
        genes.loc[i, ["role_hypothesis", "role_rule_id", "next_experiment_rule_id"]] = [role, role_rule, experiment]
        genes.loc[i, "allowed_wording"] = "The genetic domain is untestable under the adopted priority eligibility rule; RNA evidence remains separate."
        genes.loc[i, "limitation"] += " Genetics eligibility was corrected by the approved priority rule; an untestable genetic domain cannot supply a negative test."
    for frame in (genes, evidence, coverage):
        frame["analysis_release_id"] = rel
        frame["passport_id"] = rel + ":" + frame.ensembl_id.astype(str)
    genetic_mask = evidence.evidence_domain.eq("genetics") & evidence.ensembl_id.isin(moves.index)
    evidence.loc[genetic_mask, "call_state"] = "untestable"
    evidence.loc[genetic_mask, ["testability_state", "testability_reason"]] = "insufficient_shared_posterior"
    evidence.loc[genetic_mask, "assay"] = "SuSiE_COLOC"
    # Historical ABF/SuSiE reference values remain in the method companion,
    # never as the result of an untestable support call.
    for column in ("estimate", "standard_error", "ci_lower", "ci_upper", "p_value", "q_value"):
        evidence.loc[genetic_mask, column] = pd.NA
    evidence.loc[genetic_mask, "allowed_wording"] = "SuSiE-COLOC untestable under the adopted priority eligibility rule; not a genetic negative."
    evidence.loc[genetic_mask, "limitation"] += " Original display values are retained only as labeled method-specific reference evidence."
    cov_mask = coverage.evidence_domain.eq("genetics") & coverage.ensembl_id.isin(moves.index)
    coverage.loc[cov_mask, "call_state"] = "untestable"
    coverage.loc[cov_mask, ["testability_state", "testability_reason"]] = "insufficient_shared_posterior"
    coverage.loc[cov_mask, "assay"] = "SuSiE_COLOC"
    cov_original = coverage.evidence_domain.eq("genetics")
    coverage.loc[cov_original, "assay"] = coverage.loc[cov_original, "ensembl_id"].map(evidence[evidence.evidence_domain.eq("genetics")].set_index("ensembl_id").assay)
    for name, frame in (("passport_gene_index.parquet", genes), ("passport_evidence_long.parquet", evidence), ("passport_coverage_long.parquet", coverage)):
        frame.to_parquet(out / name, index=False)
    for name in ("passport_program_index.parquet", "passport_program_membership.parquet", "passport_program_context.parquet"):
        frame = pd.read_parquet(source / name)
        frame["analysis_release_id"] = rel
        frame.to_parquet(out / name, index=False)
    for name in ("passport_domain_coverage.tsv", "passport_assay_status.tsv", "passport_dataset_status.tsv", "passport_source_nodes.tsv", "passport_source_edges.tsv", "passport_rulebook.tsv", "passport_controlled_vocabularies.tsv", "passport_data_dictionary.tsv", "passport_identity_quarantine.tsv"):
        frame = tsv(source / name)
        if "analysis_release_id" in frame:
            frame["analysis_release_id"] = rel
        save_tsv(frame, out / name)
    experiments = tsv(source / "passport_next_experiment.tsv").set_index("ensembl_id")
    for ens in changed_gene_ids:
        experiments.loc[ens, "experiment_rule_id"] = "EXP_DIRECT_ASSAY_V1"
        experiments.loc[ens, "role_rule_id"] = "ROLE_UNRESOLVED_V1"
        for key, value in contract.EXPERIMENT_RULES["EXP_DIRECT_ASSAY_V1"].items():
            experiments.loc[ens, key] = value
    experiments["analysis_release_id"] = rel
    experiments["passport_id"] = rel + ":" + experiments.index
    save_tsv(experiments.reset_index(), out / "passport_next_experiment.tsv")
    save_tsv(moves.reset_index(), out / "eligibility_changes.tsv")
    save_tsv(old_genes[["ensembl_id", "symbol", "primary_evidence_class", "next_experiment_rule_id"]].rename(columns={"primary_evidence_class": "named_r8_class", "next_experiment_rule_id": "named_r8_next_experiment_rule"}), out / "named_r8_reference.tsv")

    gmap = genes.set_index("ensembl_id")
    current_gen = evidence[evidence.evidence_domain == "genetics"].set_index("ensembl_id")
    current_rna = evidence[evidence.evidence_domain == "transcriptomics"].set_index("ensembl_id")
    method = tsv(companion / "genetics_method_evidence.tsv")
    method = method.rename(columns={"primary_evidence_class": "named_r8_class", "r8_genetics_call_state": "named_r8_genetics_call_state"})
    method["primary_evidence_class"] = method.ensembl_id.map(gmap.primary_evidence_class)
    method["genetics_call_state"] = method.ensembl_id.map(current_gen.call_state)
    method["adoption_status"] = "adopted"
    save_tsv(method, out / "genetics_method_evidence.tsv")
    shutil.copy2(root / c["genetics_qualification"], out / "genetics_qualification.tsv")

    molecular_records = []
    with inputs["corrected_molecular"].open() as handle, (out / "molecular_object_sidecar.jsonl").open("x") as dest:
        for line in handle:
            record = json.loads(line)
            ens = record["ensembl_id_base"]
            record["named_r8_class"] = record["primary_evidence_class"]
            record["primary_evidence_class"] = gmap.loc[ens, "primary_evidence_class"]
            record["genetics_call"] = current_gen.loc[ens, "call_state"]
            record["legacy_next_experiment_rule_id"] = gmap.loc[ens, "next_experiment_rule_id"]
            if ens in changed_gene_ids:
                record["candidate_object"] = "untestable_nomination"
                record["open_mechanistic_question"] = "missing_decisive_assay"
            record["analysis_release_id"] = rel
            record["catalog_object_id"] = rel + ":" + ens + ":" + record["candidate_object"]
            record["candidate_not_adopted"] = False
            record["adoption_status"] = "adopted"
            record["input_state"] = "adopted_priority_eligibility_with_explicit_assay_coverage"
            record["recommended_experiment_rule_id"] = molecular.route_experiment(record)
            record["next_experiment_rule_id"] = record["recommended_experiment_rule_id"]
            record = molecular.adapt_record(record)
            assert record["rna_call"] == current_rna.loc[ens, "call_state"]
            if ens in changed_gene_ids:
                assert record["recommended_experiment_rule_id"] == "EXP_MEASURE_MISSING_ASSAY_V2"
            molecular_records.append(record)
            dest.write(json.dumps(record, sort_keys=True, allow_nan=False) + "\n")
    annotated = pd.read_parquet(inputs["annotation_r9d"])
    assert set(annotated.ensembl_id) == set(genes.ensembl_id)
    for col in ("primary_evidence_class", "role_hypothesis", "role_rule_id", "next_experiment_rule_id", "allowed_wording", "limitation"):
        annotated["named_r8_" + col] = annotated[col]
        annotated[col] = annotated.ensembl_id.map(gmap[col])
    annotated["analysis_release_id"] = rel
    annotated["passport_id"] = rel + ":" + annotated.ensembl_id
    annotated["genetics_call_state"] = annotated.ensembl_id.map(current_gen.call_state)
    annotated["genetics_testability_state"] = annotated.ensembl_id.map(current_gen.testability_state)
    annotated["named_r8_gen_untestable_reason_r9"] = annotated.gen_untestable_reason_r9
    annotated.loc[annotated.ensembl_id.isin(moves.index), "gen_untestable_reason_r9"] = "priority_shared_posterior_untestable"
    annotated["native_susie_state_t12"] = annotated.gen_susie_state_t12
    annotated.loc[annotated.ensembl_id.isin(moves.index), "gen_susie_state_t12"] = "untestable"
    annotated["class_display_r9"] = annotated.primary_evidence_class.map({"concordant": "both supported", "genetically_anchored": "genetically anchored", "established_state_associated": "established state-associated", "unresolved": "unresolved", "untested": "untested"})
    annotated["untested_missing_domain"] = ""
    annotated["untested_one_domain_supported"] = ""
    for i, row in annotated.iterrows():
        ens = row.ensembl_id
        if row.primary_evidence_class == "untested":
            missing = [domain for domain, frame in (("genetics", current_gen), ("RNA", current_rna)) if frame.loc[ens, "call_state"] == "untestable"]
            supported = [domain for domain, frame in (("genetics", current_gen), ("RNA", current_rna)) if frame.loc[ens, "call_state"] == "supported"]
            annotated.loc[i, "untested_missing_domain"] = ";".join(missing)
            annotated.loc[i, "untested_one_domain_supported"] = ";".join(supported)
    annotated.to_parquet(out / "catalog_annotated_gene_table.parquet", index=False)
    annotated.to_csv(out / "catalog_annotated_gene_table.tsv.gz", sep="\t", index=False, compression={"method": "gzip", "mtime": 0})
    annotations = out / "annotations"
    annotations.mkdir(exist_ok=True)
    sidecar_inputs = list(annotation.glob("sidecar_*.tsv*")) + [root / c["stage_replication"], root / c["anchor_disease"]]
    annotation_summary = []
    for path in sidecar_inputs:
        if path.name == "sidecar_h86_protein_coverage_by_class.tsv":
            continue  # Recompute this aggregate below from unchanged protein rows.
        frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        if "ensembl_id" in frame:
            ids = frame.ensembl_id.str.split(".").str[0]
            for col in ("primary_evidence_class", "class", "catalog_class", "role_hypothesis"):
                if col in frame:
                    if col == "class" and not set(frame[col]) <= set(c["expected_class_counts"]) | {""}:
                        # Figure 5B's 'class' is RNA/protein significance context,
                        # not the Catalog evidence class.
                        frame["catalog_primary_evidence_class"] = ids.map(gmap.primary_evidence_class)
                        continue
                    frame["named_r8_" + col] = frame[col]
                    target = "role_hypothesis" if col == "role_hypothesis" else "primary_evidence_class"
                    mapped = ids.map(gmap[target])
                    frame[col] = mapped.where(mapped.notna(), frame[col])
            if "genetic_call_state" in frame:
                frame["named_r8_genetic_call_state"] = frame.genetic_call_state
                mapped = ids.map(current_gen.call_state)
                frame.genetic_call_state = mapped.where(mapped.notna(), frame.genetic_call_state)
            for col in ("genetic_testability_state", "genetic_testability_reason"):
                if col in frame:
                    frame["named_r8_" + col] = frame[col]
                    mapped = ids.map(current_gen["testability_state" if col.endswith("state") else "testability_reason"])
                    frame[col] = mapped.where(mapped.notna(), frame[col])
            if "susie_state_t12" in frame:
                frame["native_susie_state_t12"] = frame.susie_state_t12
                frame.loc[ids.isin(moves.index), "susie_state_t12"] = "untestable"
            frame["catalog_release_id"] = rel
            state = "current classes joined by stable Ensembl ID; native effects and assay states unchanged"
        else:
            state = "historical source context; no current gene-class inference"
            frame["catalog_context_state"] = state
        target = annotations / path.name
        frame.to_csv(target, sep="\t", index=False, compression={"method": "gzip", "mtime": 0} if path.suffix == ".gz" else None)
        annotation_summary.append({"file": target.name, "source": str(path), "source_sha256": sha(path), "interpretation": state})
    protein = tsv(root / c["protein_assignment"])
    protein["named_r8_class"] = protein["class"]
    symbol_class = genes.set_index("symbol").primary_evidence_class
    assert symbol_class.index.is_unique
    protein["class"] = protein.catalog_symbol.map(symbol_class).fillna(protein["class"])
    save_tsv(protein, annotations / "protein_row_assignment_current.tsv")
    protein_summary = []
    for cls, group in annotated.groupby("primary_evidence_class"):
        assigned = protein[protein["class"] == cls]
        n = len(assigned)
        protein_summary.append({"class": cls, "n_genes": len(group), "n_quantified": n, "pct_quantified": 100 * n / len(group),
                                "n_quantified_direct": int(assigned.route.eq("direct").sum()),
                                "n_quantified_remapped": int(assigned.route.ne("direct").sum()),
                                "n_detected_in_deposit": int(group.liver_proteome_any.fillna(0).astype(float).sum())})
    save_tsv(pd.DataFrame(protein_summary), annotations / "protein_coverage_by_current_class.tsv")
    save_tsv(pd.DataFrame(annotation_summary), annotations / "annotation_sources.tsv")
    graph = out / "portal_export"
    graph.mkdir()
    shutil.copy2(source / "portal_export/passport_source_graph.json", graph / "passport_source_graph.json")
    for name in ("passport_gene_index.parquet", "passport_evidence_long.parquet", "passport_coverage_long.parquet"):
        shutil.copy2(out / name, graph / name)
    ui.build_ui(out)
    browser = out / "portal_candidate/index.html"
    page = browser.read_text().replace("MASLD Gene Catalog — local candidate review", "MASLD Gene Catalog — adopted eligibility release")
    page = page.replace("Candidate review · local only", "Adopted eligibility release · local access")
    page = page.replace("Candidate review artifact; no canonical promotion is authorized.", "Adopted eligibility release; local access only.")
    browser.write_text(page)
    ui_contract = out / "portal_candidate/ui_contract.json"
    uc = json.loads(ui_contract.read_text())
    uc.update(adoption_status="adopted", canonical_promotion_authorized=True, candidate_not_adopted=False)
    ui_contract.write_text(json.dumps(uc, indent=2) + "\n")

    # Biological invariants use original source tables, separately from new IDs.
    assert len(genes) == 28233 and genes.ensembl_id.tolist() == old_genes.ensembl_id.tolist()
    counts = genes.primary_evidence_class.value_counts().to_dict()
    assert counts == c["expected_class_counts"]
    pd.testing.assert_frame_equal(old_evidence[old_evidence.evidence_domain == "transcriptomics"].drop(columns=["analysis_release_id", "passport_id"]).reset_index(drop=True), evidence[evidence.evidence_domain == "transcriptomics"].drop(columns=["analysis_release_id", "passport_id"]).reset_index(drop=True))
    old_positive = set(old_genetics.index[old_genetics.call_state == "supported"])
    assert old_positive == set(current_gen.index[current_gen.call_state == "supported"]) and len(old_positive) == 436
    assert int(genetic_mask.sum()) == 188 and int(gene_mask.sum()) == 180
    pd.testing.assert_frame_equal(pd.read_parquet(inputs["source_membership"]).drop(columns="analysis_release_id"), pd.read_parquet(out / "passport_program_membership.parquet").drop(columns="analysis_release_id"))
    assert all(sha(path) == before_hash[k] for k, path in inputs.items())
    matrix = pd.crosstab(current_gen.call_state, current_rna.call_state).astype(int)
    matrix.index.name = "genetics_call"; matrix.columns.name = "RNA_call"
    matrix.to_csv(out / "call_matrix.tsv", sep="\t")
    metadata = {"release_id": rel, "adoption_status": "adopted", "candidate_not_adopted": False, "authorization": c["authorization"],
                "gene_records": len(genes), "evidence_rows": len(evidence), "programs": 117, "memberships": 6633,
                "classes": counts, "genetics_state_changes": 188, "class_changes": 180,
                "jointly_classified": int((genes.primary_evidence_class != "untested").sum()),
                "genetics_call_counts": current_gen.call_state.value_counts().to_dict(), "RNA_call_counts": current_rna.call_state.value_counts().to_dict(),
                "call_matrix": {str(k): {str(q): int(v) for q, v in vals.items()} for k, vals in matrix.to_dict().items()},
                "gene_route_counts": genes.next_experiment_rule_id.value_counts().to_dict(),
                "molecular_route_counts": dict(Counter(r["recommended_experiment_rule_id"] for r in molecular_records)),
                "RNA_effects_calls_directions_unchanged": True, "positive_genetic_membership_unchanged": True,
                "program_weights_memberships_unchanged": True, "named_inputs_unchanged": True,
                "source_inputs": {k: {"path": str(path), "sha256": before_hash[k]} for k, path in inputs.items()},
                "limitation": "Eligibility correction only; no new genetic or RNA effect, no LD/prior sensitivity membership changes, no external publication."}
    (out / "adoption_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    shutil.copy2(a.config, out / "build_config.json")
    pd.DataFrame([{"relative_path": str(path.relative_to(out)), "sha256": sha(path), "bytes": path.stat().st_size} for path in sorted(out.rglob("*")) if path.is_file()]).to_csv(out / "release_files.tsv", sep="\t", index=False)
    links = [("Browse all 28,233 genes", "portal_candidate/index.html"), ("Gene index (Parquet)", "passport_gene_index.parquet"), ("Assay evidence (Parquet)", "passport_evidence_long.parquet"), ("Enriched annotations and RNA uncertainty", "catalog_annotated_gene_table.tsv.gz"), ("Separate SuSiE and ABF posteriors, studies and ancestry", "genetics_method_evidence.tsv"), ("Molecular-object records and next experiments", "molecular_object_sidecar.jsonl"), ("Eligibility changes", "eligibility_changes.tsv"), ("Gene-level next experiments", "passport_next_experiment.tsv"), ("Counts, unchanged measurements and sources", "adoption_metadata.json")]
    (out / "index.html").write_text('<meta charset="utf-8"><title>Adopted MASLD Gene Catalog</title><style>body{font:16px sans-serif;max-width:1000px;margin:2rem auto;line-height:1.5}li{margin:.7rem 0}</style><h1>MASLD Gene Catalog</h1><p>Adopted priority-eligibility release, 2026-11-06. 28,233 genes: 13,731 untested; 12,958 unresolved; 1,118 established-state associated; 400 genetically anchored; 26 supported in both domains. Gene classes preserve separate RNA and genetic support; directions were not compared for the both-supported class.</p><p>188 genetics states and 180 classes are corrected under the approved priority rule. All RNA effects and calls, all positive genetic memberships and the 117 fixed program weights are unchanged. Native SuSiE/ABF reference posteriors are displayed separately. Untestable evidence is not a negative result. The annotated download retains RNA standard errors, nominal P values, BH q values and the original r8 fields.</p><ul>' + ''.join('<li><a href="'+html.escape(url)+'">'+html.escape(label)+'</a></li>' for label,url in links) + '</ul><p>This is local file access; no external portal has been published. Historical r8/r9 outputs remain preserved.</p>')
    print(json.dumps(metadata))


if __name__ == "__main__":
    main()
