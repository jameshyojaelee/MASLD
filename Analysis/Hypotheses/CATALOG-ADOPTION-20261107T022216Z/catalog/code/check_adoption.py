"""Independent whole-Catalog comparison against pre-adoption tables and approved IDs."""
from collections import Counter
import hashlib
import json
from pathlib import Path
import shutil

import pandas as pd

root = Path(__file__).resolve().parents[5]
work = root / "Analysis/Hypotheses/CATALOG-ADOPTION-20261107T022216Z/catalog"
config_path = root / "scripts/portal/catalog_current.json"
config = json.loads(config_path.read_text())
out = root / config["output"]
old = root / config["source_bundle"]
read = lambda path: pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
g0, g = [pd.read_parquet(p / "passport_gene_index.parquet").set_index("ensembl_id") for p in (old, out)]
e0, e = [pd.read_parquet(p / "passport_evidence_long.parquet").set_index(["ensembl_id", "evidence_domain"]) for p in (old, out)]
cv = pd.read_parquet(out / "passport_coverage_long.parquet").set_index(["ensembl_id", "evidence_domain"])
proposal = read(root / config["eligibility_proposals"])
approved = set(proposal.loc[proposal.priority_rule_eligibility_transition == "True", "ensembl"])
class_approved = set(proposal.loc[proposal.priority_rule_class_would_change == "True", "ensembl"])
assert len(approved) == 188 and len(class_approved) == 180
assert g.index.tolist() == g0.index.tolist() and len(g) == 28233 and g.index.is_unique
pd.testing.assert_series_equal(g.symbol, g0.symbol)
assert set(g.index[g.primary_evidence_class != g0.primary_evidence_class]) == class_approved
assert g.primary_evidence_class.value_counts().to_dict() == config["expected_class_counts"]
changed_genetics = {ens for ens in g.index if e.loc[(ens, "genetics"), "call_state"] != e0.loc[(ens, "genetics"), "call_state"]}
assert changed_genetics == approved
for ens in approved:
    x, coverage = e.loc[(ens, "genetics")], cv.loc[(ens, "genetics")]
    assert x.call_state == coverage.call_state == "untestable"
    assert x.testability_state == coverage.testability_state == "insufficient_shared_posterior"
    assert pd.isna(x.estimate) and pd.isna(x.q_value) and pd.isna(x.p_value)
for domain in ("transcriptomics",):
    before, after = e0.xs(domain, level="evidence_domain"), e.xs(domain, level="evidence_domain")
    cols = [x for x in before if x not in {"analysis_release_id", "passport_id"}]
    pd.testing.assert_frame_equal(before[cols], after[cols])
assert set(e0.xs("genetics", level="evidence_domain").query("call_state == 'supported'").index) == set(e.xs("genetics", level="evidence_domain").query("call_state == 'supported'").index)
assert int(e.xs("genetics", level="evidence_domain").call_state.eq("supported").sum()) == 436
pd.testing.assert_frame_equal(pd.read_parquet(old / "passport_program_membership.parquet").drop(columns="analysis_release_id"), pd.read_parquet(out / "passport_program_membership.parquet").drop(columns="analysis_release_id"))
pd.testing.assert_frame_equal(pd.read_parquet(old / "passport_program_index.parquet").drop(columns="analysis_release_id"), pd.read_parquet(out / "passport_program_index.parquet").drop(columns="analysis_release_id"))
assert len(e) == len(cv) == 56466
assert e.call_state.equals(cv.call_state.reindex(e.index))
experiments = read(out / "passport_next_experiment.tsv").set_index("ensembl_id")
assert set(experiments.index) == set(g.index)
assert experiments.loc[list(class_approved), "experiment_rule_id"].eq("EXP_DIRECT_ASSAY_V1").all()
assert g.loc[list(class_approved), "next_experiment_rule_id"].eq("EXP_DIRECT_ASSAY_V1").all()
domain = read(out / "passport_domain_coverage.tsv").set_index("evidence_domain")
assert int(domain.loc["genetics", "n_gene_evidence_rows"]) == int(domain.loc["transcriptomics", "n_gene_evidence_rows"]) == 28233
assert int(domain.loc["cell_context", "n_programs"]) == 117
assert int(domain.loc["cell_context", "n_program_membership_rows"]) == 6633
records = []
with (out / "molecular_object_sidecar.jsonl").open() as handle:
    records = [json.loads(line) for line in handle]
assert len(records) == 28233 and len({r["ensembl_id_base"] for r in records}) == 28233
for r in records:
    ens = r["ensembl_id_base"]
    assert r["primary_evidence_class"] == g.loc[ens, "primary_evidence_class"]
    assert r["genetics_call"] == e.loc[(ens, "genetics"), "call_state"]
    assert r["rna_call"] == e.loc[(ens, "transcriptomics"), "call_state"]
    assert r["candidate_not_adopted"] is False and r["adoption_status"] == "adopted"
    if r["rna_call"] == "untestable":
        assert json.loads(r["assay_applicability"])["bulk_rna"] in {"untestable", "not_applicable"}
    if ens in class_approved:
        assert r["recommended_experiment_rule_id"] == "EXP_MEASURE_MISSING_ASSAY_V2"
annotation0 = pd.read_parquet(root / config["annotation_bundle"] / "catalog_r9d_gene_table.parquet").set_index("ensembl_id")
annotation = pd.read_parquet(out / "catalog_annotated_gene_table.parquet").set_index("ensembl_id")
rna_cols = [c for c in annotation0 if c.startswith("rna_") or c.startswith("fig3e_") or c.startswith("h51_") or c.startswith("coding_cs_")]
pd.testing.assert_frame_equal(annotation0[rna_cols], annotation[rna_cols])
pd.testing.assert_series_equal(annotation.primary_evidence_class, g.primary_evidence_class)
for ens in approved:
    assert annotation.loc[ens, "genetics_call_state"] == "untestable"
    assert annotation.loc[ens, "gen_susie_state_t12"] == "untestable"
protein_source = read(root / config["protein_assignment"])
protein = read(out / "annotations/protein_row_assignment_current.tsv")
assert len(protein) == len(protein_source) == 6071
pd.testing.assert_frame_equal(protein_source.drop(columns="class"), protein.drop(columns=["class", "named_r8_class"]))
protein_expected = protein_source.catalog_symbol.map(g.set_index("symbol").primary_evidence_class)
assert (protein.loc[protein_expected.notna(), "class"] == protein_expected[protein_expected.notna()]).all()
coverage = read(out / "annotations/protein_coverage_by_current_class.tsv").set_index("class")
assert coverage.n_quantified.astype(int).to_dict() == {"concordant": 11, "established_state_associated": 170, "genetically_anchored": 203, "unresolved": 5171, "untested": 487}
# The Figure 5B significance-context class is a different measurement.
p0 = read(root / config["annotation_bundle"] / "sidecar_fig5b_liver_protein.tsv")
p1 = read(out / "annotations/sidecar_fig5b_liver_protein.tsv")
pd.testing.assert_frame_equal(p0, p1[p0.columns])
for path in (out / "annotations").glob("*"):
    if path.name in {"annotation_sources.tsv", "protein_row_assignment_current.tsv", "protein_coverage_by_current_class.tsv"}:
        continue
    source = root / config["annotation_bundle"] / path.name
    if path.name == Path(config["stage_replication"]).name: source = root / config["stage_replication"]
    if path.name == Path(config["anchor_disease"]).name: source = root / config["anchor_disease"]
    if not source.exists(): continue
    before, after = read(source), read(path)
    # Preserve every native effect, SE, p/q, direction and biological assay state.
    preserved = [k for k in before if any(token in k.lower() for token in ("logfc", "log2fc", "coef", "slope", "padj", "pvalue", "direction", "_se", "_ci", "protein_log", "protein_p"))]
    if preserved: pd.testing.assert_frame_equal(before[preserved], after[preserved])
metadata = json.loads((out / "adoption_metadata.json").read_text())
assert metadata["adoption_status"] == "adopted" and metadata["candidate_not_adopted"] is False
assert metadata["jointly_classified"] == 14502
for value in metadata["source_inputs"].values():
    assert sha(Path(value["path"])) == value["sha256"]
report = {"status": "PASS", "records": 28233, "evidence_rows": 56466, "genetics_state_changes": 188, "class_changes": 180,
          "joint_eligible": 14502, "all_RNA_effects_calls_directions_SE_p_q_unchanged": True,
          "positive_genetics_unchanged": 436, "programs_unchanged": 117, "membership_weights_unchanged": 6633,
          "all_180_gene_and_molecular_routes_corrected": True, "no_unavailable_genetic_statistic_reported_as_a_test": True,
          "complete_protein_source_rows_preserved": 6071, "protein_quantified_by_class": coverage.n_quantified.astype(int).to_dict(),
          "annotation_native_values_preserved": True, "domain_coverage_totals_correct": True,
          "molecular_RNA_untestability_correct": True, "named_inputs_unchanged": True}
(work / "checks/independent_adoption_checks.json").write_text(json.dumps(report, indent=2) + "\n")
(out / "independent_checks.json").write_text(json.dumps(report, indent=2) + "\n")
producers = out / "producers"; producers.mkdir(exist_ok=True)
for path in (root / "scripts/portal/build_adopted_eligibility_catalog.py", config_path, Path(__file__)):
    shutil.copy2(path, producers / path.name)
manifest = [{"relative_path": str(path.relative_to(out)), "sha256": sha(path), "bytes": path.stat().st_size} for path in sorted(out.rglob("*")) if path.is_file() and path.name != "release_files.tsv"]
pd.DataFrame(manifest).to_csv(out / "release_files.tsv", sep="\t", index=False)
print(json.dumps(report))
