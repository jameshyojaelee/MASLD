#!/usr/bin/env python3
"""Independent checks of the approved Catalog change, against the old full tables."""
import json
from collections import Counter
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
RUN = Path(__file__).resolve().parent
OLD = ROOT / "figures/candidates/catalog-r8-completion-20261014T141645Z_22236527/catalog_local"
NEW = ROOT / "RNA-seq/results/evidence_passports/releases/gene-catalog-eligibility-20261106"
PROPOSAL = ROOT / "Analysis/Hypotheses/REVIEW-FIX-20261107T020001Z/catalog_models/outputs/eligibility_candidate/catalog_priority_eligibility_candidate.tsv"


def read(bundle, name):
    return pd.read_parquet(bundle / name)


def main():
    old = read(OLD, "passport_gene_index.parquet").set_index("ensembl_id").sort_index()
    new = read(NEW, "passport_gene_index.parquet").set_index("ensembl_id").sort_index()
    proposal = pd.read_csv(PROPOSAL, sep="\t", dtype=str).set_index("ensembl_id").sort_index()
    assert old.index.tolist() == new.index.tolist() == proposal.index.tolist()
    proposal.index = proposal.index.astype(new.index.dtype)
    assert new.index.is_unique and len(new) == 28233
    pd.testing.assert_series_equal(new.primary_evidence_class, proposal.candidate_priority_class, check_names=False, check_dtype=False)
    changed = old.primary_evidence_class.ne(new.primary_evidence_class)
    assert int(changed.sum()) == 180
    assert old.loc[changed, "primary_evidence_class"].value_counts().to_dict() == {"unresolved": 169, "established_state_associated": 11}
    oe = read(OLD, "passport_evidence_long.parquet").set_index(["ensembl_id", "evidence_domain"]).sort_index()
    ne = read(NEW, "passport_evidence_long.parquet").set_index(["ensembl_id", "evidence_domain"]).sort_index()
    assert len(ne) == 56466 and oe.index.equals(ne.index) and ne.index.is_unique
    orna = oe.xs("transcriptomics", level="evidence_domain")
    nrna = ne.xs("transcriptomics", level="evidence_domain")
    rna_columns = [c for c in orna if c not in {"analysis_release_id", "passport_id"}]
    pd.testing.assert_frame_equal(orna[rna_columns], nrna[rna_columns], check_dtype=False)
    og = oe.xs("genetics", level="evidence_domain")
    ng = ne.xs("genetics", level="evidence_domain")
    genetic_changes = og.call_state.ne(ng.call_state)
    assert int(genetic_changes.sum()) == 188
    assert set(ng.loc[genetic_changes, "call_state"]) == {"untestable"}
    assert set(og.loc[genetic_changes, "call_state"]) == {"indeterminate"}
    assert set(ng.index[ng.call_state.eq("supported")]) == set(og.index[og.call_state.eq("supported")])
    assert int(ng.call_state.eq("supported").sum()) == 436
    classes = []
    for ens in new.index:
        g, r = ng.loc[ens, "call_state"], nrna.loc[ens, "call_state"]
        if "untestable" in (g, r):
            cls = "untested"
        elif g == "supported" and r == "supported":
            cls = "concordant"
        elif g == "supported":
            cls = "genetically_anchored"
        elif r == "supported":
            cls = "established_state_associated"
        else:
            cls = "unresolved"
        assert new.loc[ens, "primary_evidence_class"] == cls
        classes.append(cls)
    assert sum(c != "untested" for c in classes) == 14502
    for filename in ("passport_program_index.parquet", "passport_program_membership.parquet", "passport_program_context.parquet"):
        a, b = read(OLD, filename), read(NEW, filename)
        cols = [c for c in a if c != "analysis_release_id"]
        pd.testing.assert_frame_equal(a[cols], b[cols])
    routes = pd.read_csv(NEW / "passport_next_experiment.tsv", sep="\t").set_index("ensembl_id")
    assert (routes.loc[new.index[changed], "experiment_rule_id"] == "EXP_DIRECT_ASSAY_V1").all()
    molecular_counts = Counter()
    records = 0
    with (NEW / "molecular_object_sidecar.jsonl").open() as handle:
        for line in handle:
            row = json.loads(line)
            ens = row["ensembl_id_base"]
            assert row["primary_evidence_class"] == new.loc[ens, "primary_evidence_class"]
            assert row["rna_call"] == nrna.loc[ens, "call_state"]
            assert row["genetics_call"] == ng.loc[ens, "call_state"]
            assert not row["candidate_not_adopted"]
            if changed.loc[ens]:
                assert row["recommended_experiment_rule_id"] == "EXP_MEASURE_MISSING_ASSAY_V2"
            molecular_counts[row["recommended_experiment_rule_id"]] += 1
            records += 1
    assert records == 28233
    annotated = read(NEW, "catalog_annotated_gene_table.parquet").set_index("ensembl_id").sort_index()
    pd.testing.assert_series_equal(annotated.primary_evidence_class, new.primary_evidence_class)
    result = {"status": "PASS", "genes": len(new), "evidence_rows": len(ne), "genetics_state_changes": int(genetic_changes.sum()), "class_changes": int(changed.sum()), "class_counts": dict(Counter(classes)), "all_RNA_fields_unchanged_except_release_identifiers": True, "genetic_positive_membership_unchanged": True, "program_definitions_weights_memberships_unchanged": True, "approved_candidate_classes_match": True, "classes_rederived_from_domain_states": True, "changed_routes_consistent": True, "molecular_route_counts": dict(molecular_counts)}
    (RUN / "independent_adoption_checks.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
