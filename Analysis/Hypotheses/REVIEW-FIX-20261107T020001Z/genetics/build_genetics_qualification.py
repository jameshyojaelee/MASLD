#!/usr/bin/env python3
"""Expose saved genetic uncertainty without modifying named Catalog classes.

Standard-library-only, small-table reconstruction. No GWAS fitting or biological
inference is run. Output must be a new directory. Alternative classes are proposals.
"""
import argparse
import csv
import hashlib
import json
import math
import platform
from collections import Counter
from pathlib import Path


def read_table(path, delimiter="\t"):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def finite(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def changed_prior(pp4, ratio):
    """Reweight H4 only; preserve a missing result as missing."""
    if pp4 is None:
        return ""
    if not 0 <= pp4 <= 1 or ratio <= 0:
        raise ValueError("invalid posterior or prior ratio")
    return ratio * pp4 / (ratio * pp4 + (1 - pp4))


def write_table(path, rows):
    if not rows:
        raise ValueError(f"empty output: {path}")
    with path.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def keyed(rows, key):
    result = {row[key]: row for row in rows}
    if len(result) != len(rows):
        raise ValueError(f"duplicate {key}")
    return result


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--catalog-identities", type=Path,
                        help="Optional checked Catalog identity TSV; joins only on stable ID")
    args = parser.parse_args()
    root = args.root.resolve()
    review = root / "Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses"
    an1 = review / "AN1/gene_level_20261021T165252Z"
    r2 = root / "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z"
    sources = {
        "native_tier12": r2 / "gene_level_coloc_tier12.csv",
        "eligibility_states": an1 / "gene_states.tsv",
        "priority_catalog_moves": an1 / "catalog_moves_rule_a.tsv",
        "any_eligible_catalog_moves": an1 / "catalog_moves_rule_b.tsv",
        "eligibility_summary": an1 / "results.json",
        "ld_diagnostics": review / "AN5/pph4_1_diagnostic_summary.tsv",
        "ld_sole_support_roster": review / "AN5/pph4_1_rows_positive_t12.tsv",
        "saved_prior_summary": review / "AN5/p12/p12_summary.json",
        "saved_catalog_prior_support": review / "AN5/p12/p12_per_gene_supported_or_gained.tsv",
        "replication_gene_records": root / "Analysis/Hypotheses/RV9-20261102T151012Z/out-20261102T151110Z/gene_states.tsv",
        "builder": Path(__file__).resolve(),
    }
    if args.catalog_identities:
        sources["catalog_identities"] = args.catalog_identities.resolve()
    before = {name: digest(path) for name, path in sources.items()}
    native = keyed(read_table(sources["native_tier12"], ","), "ensembl")
    states = keyed(read_table(sources["eligibility_states"]), "ensembl")
    moves = {rule: keyed(read_table(sources[name]), "ensembl") for rule, name in
             [("a", "priority_catalog_moves"), ("b", "any_eligible_catalog_moves")]}
    assert native.keys() == states.keys(), "native and saved sensitivity identities differ"
    names = {key: row["gene"] for key, row in
             keyed(read_table(sources["replication_gene_records"]), "ensembl").items()}
    catalog = {}
    if args.catalog_identities:
        catalog = keyed(read_table(sources["catalog_identities"]), "ensembl_id")
        assert len(catalog) == 28233
        names.update({key: row["symbol"] for key, row in catalog.items()})
    for mapping in moves.values():
        names.update({key: row["symbol"] for key, row in mapping.items()})

    diagnostics = read_table(sources["ld_diagnostics"])
    roster = read_table(sources["ld_sole_support_roster"])
    diagnostic_ids = {}
    for row in roster:
        previous = diagnostic_ids.setdefault(row["gene"], row["ensembl"])
        assert previous == row["ensembl"], "ambiguous diagnostic gene identity"
        names[row["ensembl"]] = row["gene"]
    ld_by_gene = {}
    for row in diagnostics:
        if row["role"] == "pph4_1":
            ld_by_gene.setdefault(row["gene"], []).append(row)
    assert len(ld_by_gene) == 14 and sum(map(len, ld_by_gene.values())) == 18
    assert set(ld_by_gene) == set(diagnostic_ids)
    ld_rows = []
    for symbol, fits in sorted(ld_by_gene.items()):
        values = [finite(row["pp4_drop_dentist"]) for row in fits]
        finite_values = [value for value in values if value is not None]
        max_pp4 = max(finite_values) if finite_values else None
        ld_rows.append({
            "ensembl": diagnostic_ids[symbol], "symbol": symbol,
            "n_near_certain_sole_support_fits": len(fits),
            "n_anchor_abs_z_residual_gt_5_45": sum(abs(float(row["z_std_diff_at_best_hit1"])) > 5.45 for row in fits),
            "max_susie_pp4_after_posthoc_variant_removal": "" if max_pp4 is None else max_pp4,
            "retains_any_tested_fit_above_0_5": bool(max_pp4 is not None and max_pp4 > 0.5),
            "posthoc_sensitivity_state": "no_credible_set_result" if not finite_values else
                "retains_support" if max_pp4 > 0.5 else "support_below_threshold",
            "interpretation": "LD-sensitive candidate; outlier removal is a post hoc sensitivity, not proof of a false association",
            "named_class_changed": False,
        })
    ld_map = keyed(ld_rows, "ensembl")

    gene_rows = []
    for ensembl, row in sorted(native.items()):
        state = states[ensembl]
        pp4 = finite(row["coloc_best_susie_pp4"])
        abf = finite(row["coloc_best_abf_pp4"])
        assert pp4 == finite(state["coloc_best_susie_pp4_base"])
        assert abf == finite(state["coloc_best_abf_pp4_base"])
        assert row["susie_state_t12"] == state["state_a_base"]
        low, high = changed_prior(pp4, 0.2), changed_prior(pp4, 2.0)
        positive = pp4 is not None and pp4 > 0.5
        diagnostic = ld_map.get(ensembl)
        gene_rows.append({
            "ensembl": ensembl, "symbol": names.get(ensembl, ""),
            "named_catalog_class": catalog.get(ensembl, {}).get("primary_evidence_class", "not_joined_to_catalog_identity"),
            "native_susie_pp4": "" if pp4 is None else pp4,
            "native_abf_pp4": "" if abf is None else abf,
            "native_susie_driver": row["susie_driving_gwas"],
            "native_abf_driver": row["abf_driving_gwas"],
            "native_susie_state": row["susie_state_t12"],
            "priority_state_after_saved_eligibility_check": state["state_a_sens"],
            "any_eligible_state_before_check": state["state_b_base"],
            "any_eligible_state_after_saved_eligibility_check": state["state_b_sens"],
            "susie_pp4_p12_1e_6": low, "susie_pp4_p12_1e_5": high,
            "support_lost_at_p12_1e_6": bool(positive and low <= 0.5),
            "support_retained_at_p12_1e_6": bool(positive and low > 0.5),
            "support_gained_at_p12_1e_5": bool(pp4 is not None and pp4 <= 0.5 and high > 0.5),
            "ld_diagnostic_scope": "selected_near_certain_sole_support_fit" if diagnostic else "not_screened_by_this_diagnostic",
            "ld_posthoc_sensitivity_state": diagnostic["posthoc_sensitivity_state"] if diagnostic else "not_assessed",
            "named_class_changed": False,
        })

    proposals = []
    for ensembl in sorted(set(moves["a"]) | set(moves["b"])):
        selected = moves["a"].get(ensembl) or moves["b"][ensembl]
        current_class = selected["primary_evidence_class"]
        for rule in "ab":
            if ensembl in moves[rule]:
                assert moves[rule][ensembl]["primary_evidence_class"] == current_class
        proposals.append({
            "ensembl": ensembl, "symbol": selected["symbol"],
            "named_catalog_class": current_class,
            "named_catalog_genetic_evidence_state": selected["r8_genetics_call_state"],
            "priority_rule_eligibility_transition": ensembl in moves["a"],
            "proposed_class_priority_rule": "untested" if ensembl in moves["a"] else current_class,
            "priority_rule_class_would_change": ensembl in moves["a"] and current_class != "untested",
            "any_eligible_rule_eligibility_transition": ensembl in moves["b"],
            "proposed_class_any_eligible_rule": "untested" if ensembl in moves["b"] else current_class,
            "any_eligible_rule_class_would_change": ensembl in moves["b"] and current_class != "untested",
            "status": "proposed sensitivity classes only; no named Catalog record modified",
        })

    positive = {key for key, row in native.items() if (finite(row["coloc_best_susie_pp4"]) or 0) > 0.5}
    abf_positive = {key for key, row in native.items() if (finite(row["coloc_best_abf_pp4"]) or 0) > 0.5}
    assert len(native) == 18976 and len(positive) == 436 and len(abf_positive) == 822
    assert len(positive | abf_positive) == 996 and len(abf_positive - positive) == 560
    for rule in "ab":
        assert {key for key, row in states.items() if row[f"state_{rule}_sens"] == "positive"} == positive
    counts = {
        "native_tier12_genes": len(native), "susie_positive_genes": len(positive),
        "abf_positive_genes": len(abf_positive), "abf_only_genes": len(abf_positive - positive),
        "susie_or_abf_genes": len(positive | abf_positive),
        "native_susie_states": dict(Counter(row["susie_state_t12"] for row in native.values())),
        "low_prior_losses": sum(row["support_lost_at_p12_1e_6"] for row in gene_rows),
        "low_prior_retained": sum(row["support_retained_at_p12_1e_6"] for row in gene_rows),
        "high_prior_gains_native_universe": sum(row["support_gained_at_p12_1e_5"] for row in gene_rows),
        "selected_ld_diagnostic_genes": len(ld_rows),
        "selected_ld_diagnostic_fits": sum(row["n_near_certain_sole_support_fits"] for row in ld_rows),
        "ld_inconsistent_best_anchors": sum(row["n_anchor_abs_z_residual_gt_5_45"] for row in ld_rows),
        "selected_ld_genes_losing_all_tested_support_in_posthoc_sensitivity": sum(not row["retains_any_tested_fit_above_0_5"] for row in ld_rows),
        "priority_eligibility_transitions": len(moves["a"]),
        "priority_proposed_class_changes": sum(row["priority_rule_class_would_change"] for row in proposals),
        "any_eligible_eligibility_transitions": len(moves["b"]),
        "any_eligible_proposed_class_changes": sum(row["any_eligible_rule_class_would_change"] for row in proposals),
        "named_classes_modified": 0,
    }
    expected = json.loads(sources["saved_prior_summary"].read_text())
    catalog_prior = read_table(sources["saved_catalog_prior_support"])
    counts["high_prior_gains_saved_catalog_universe"] = sum(row["gain_at_p12_1e-5"] == "True" for row in catalog_prior)
    if catalog:
        gains_in_catalog = [row for row in gene_rows if row["support_gained_at_p12_1e_5"] and row["ensembl"] in catalog]
        assert len(gains_in_catalog) == counts["high_prior_gains_saved_catalog_universe"]
        counts["high_prior_gains_native_ids_absent_from_catalog"] = [row["ensembl"] for row in gene_rows
            if row["support_gained_at_p12_1e_5"] and row["ensembl"] not in catalog]
    assert counts["low_prior_losses"] == expected["lose_at_1e-6"] == 122
    assert counts["low_prior_retained"] == 314
    assert counts["selected_ld_genes_losing_all_tested_support_in_posthoc_sensitivity"] == 11
    assert counts["priority_proposed_class_changes"] == 180
    assert counts["any_eligible_proposed_class_changes"] == 39
    assert counts["high_prior_gains_saved_catalog_universe"] == expected["gain_at_1e-5"] == 29
    # Native stable IDs and annotation-adjudicated Catalog identities are distinct.
    assert counts["high_prior_gains_native_universe"] == 30
    assert changed_prior(0, 0.2) == 0 and changed_prior(1, 0.2) == 1
    assert math.isclose(changed_prior(5 / 6, 0.2), 0.5)
    assert changed_prior(None, 0.2) == ""

    args.out.mkdir(parents=True, exist_ok=False)
    write_table(args.out / "genetics_qualification.tsv", gene_rows)
    write_table(args.out / "catalog_eligibility_proposals.tsv", proposals)
    write_table(args.out / "ld_sensitive_genes.tsv", ld_rows)
    write_table(args.out / "saved_ld_fit_diagnostics.tsv", diagnostics)
    after = {name: digest(path) for name, path in sources.items()}
    assert before == after, "an input changed while building the companion"
    report = {
        "status": "validated diagnostic companion; alternative classes are proposals, not adopted",
        "counts": counts,
        "rules": {
            "priority": "positive > any untestable study > tested below threshold > no SuSiE result",
            "any_eligible": "positive > any eligible finite study > untestable > no SuSiE result",
            "prior": "same p1,p2; reweight H4 by p12 ratio; thresholds strictly above 0.5",
            "prior_universes": "30 gains in native tier1/2 gene table; 29 gains in saved annotation-adjudicated Catalog; these denominators are not interchangeable",
            "LD": "18 selected sole-support fits only; not a screen of all genetic candidates; post hoc variant removal cannot prove false support",
            "class_proposals": "only source-proven tested-to-untestable transitions; proposed untested class preserves support in other measured domains",
        },
        "checks": {"unique_native_ids": True, "saved_native_values_equal": True,
                   "positive_membership_unchanged_under_both_eligibility_rules": True,
                   "named_input_hashes_unchanged": True, "prior_boundary_and_missingness_checks": True},
        "environment": {"python": platform.python_version(), "platform": platform.platform()},
        "inputs": {name: {"path": str(sources[name]), "sha256": checksum} for name, checksum in before.items()},
        "outputs": {path.name: digest(path) for path in sorted(args.out.glob("*.tsv"))},
    }
    (args.out / "checks.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(counts, indent=2))


if __name__ == "__main__":
    main()
