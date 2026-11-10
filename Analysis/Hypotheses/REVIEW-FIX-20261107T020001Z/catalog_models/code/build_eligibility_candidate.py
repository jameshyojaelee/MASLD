#!/usr/bin/env python3
"""Make a separate, explicitly unadopted priority-eligibility Catalog view."""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import html
import json
import os
from pathlib import Path

import pandas as pd


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--project-root", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    root, out = a.project_root.resolve(), a.output.resolve()
    task = root / "Analysis/Hypotheses/REVIEW-FIX-20261107T020001Z"
    current = task / "catalog_models/outputs/corrected_companion_v2"
    proposal_path = task / "genetics/checked_companion_v2/catalog_eligibility_proposals.tsv"
    qualification_path = task / "genetics/checked_companion_v2/genetics_qualification.tsv"
    identities_path = current / "catalog_identity_current_classes.tsv"
    evidence_path = current / "passport_evidence_long_method_labeled.parquet"
    inputs = [proposal_path, qualification_path, identities_path, evidence_path]
    input_hashes = {str(path): digest(path) for path in inputs}
    out.mkdir(parents=True, exist_ok=False)
    identities = read(identities_path)
    proposals = {r["ensembl"]: r for r in read(proposal_path)}
    evidence = pd.read_parquet(evidence_path)
    rna = {str(r.ensembl_id): r._asdict() for r in evidence[evidence.evidence_domain == "transcriptomics"].itertuples(index=False)}
    assert len(identities) == len(rna) == 28233
    rows, changes = [], []
    for identity in identities:
        ens = identity["ensembl_id"]
        proposal = proposals.get(ens, {})
        state_changes = proposal.get("priority_rule_eligibility_transition", "False") == "True"
        old_class = identity["primary_evidence_class"]
        new_class = proposal["proposed_class_priority_rule"] if state_changes else old_class
        old_genetic = identity["r8_genetics_call_state"]
        new_genetic = "untestable" if state_changes else old_genetic
        if state_changes:
            assert proposal["named_catalog_class"] == old_class
            assert proposal["named_catalog_genetic_evidence_state"] == old_genetic
            assert old_genetic == "indeterminate" and new_class == "untested"
        native_rna = rna[ens]
        value = lambda key: native_rna[key] if pd.notna(native_rna[key]) else None
        row = {"ensembl_id": ens, "symbol": identity["symbol"],
               "named_r8_class": old_class, "candidate_priority_class": new_class,
               "named_r8_genetics_call": old_genetic, "candidate_priority_genetics_call": new_genetic,
               "priority_eligibility_transition": state_changes,
               "priority_class_changes": new_class != old_class,
               "rna_call_unchanged": value("call_state"), "rna_direction_unchanged": value("direction"),
               "rna_log2fc_unchanged": value("estimate"), "rna_bh_q_unchanged": value("q_value"),
               "rna_testability_unchanged": value("testability_state"),
               "rna_dataset_unchanged": value("dataset_id"), "rna_units_unchanged": value("effect_unit"),
               "candidate_not_adopted": True,
               "rule": "positive > any untestable study > tested below threshold > no SuSiE result",
               "limitation": "Priority eligibility sensitivity only; named r8 and all RNA evidence are preserved. LD and prior diagnostics remain separate annotations."}
        rows.append(row)
        if state_changes:
            changes.append(row)
    counts = dict(Counter(r["candidate_priority_class"] for r in rows))
    expected = {"untested": 13731, "unresolved": 12958, "established_state_associated": 1118, "genetically_anchored": 400, "concordant": 26}
    assert counts == expected, counts
    assert len(changes) == 188 and sum(r["priority_class_changes"] for r in rows) == 180
    assert len({r["ensembl_id"] for r in rows}) == 28233
    for name, records in (("catalog_priority_eligibility_candidate.tsv", rows), ("eligibility_changes_only.tsv", changes)):
        with (out / name).open("x") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader(); writer.writerows(records)
    # Re-read all exported RNA fields against the source, including missingness.
    reread = read(out / "catalog_priority_eligibility_candidate.tsv")
    for row in reread:
        original = rna[row["ensembl_id"]]
        for key, source_key in (("rna_call_unchanged", "call_state"), ("rna_direction_unchanged", "direction"), ("rna_testability_unchanged", "testability_state")):
            assert row[key] == original[source_key]
        for key, source_key in (("rna_log2fc_unchanged", "estimate"), ("rna_bh_q_unchanged", "q_value")):
            if pd.isna(original[source_key]):
                assert row[key] == ""
            else:
                assert float(row[key]) == float(original[source_key])
    assert all(digest(path) == input_hashes[str(path)] for path in inputs)
    checks = {"status": "PASS", "candidate_not_adopted": True, "records": len(rows),
              "genetics_state_changes": len(changes), "class_changes": sum(r["priority_class_changes"] for r in rows),
              "candidate_class_counts": counts, "named_classes_retained_in_separate_field": True,
              "all_rna_calls_directions_values_q_and_missingness_exactly_unchanged": True,
              "named_inputs_unchanged": True, "inputs": input_hashes}
    (out / "checks.json").write_text(json.dumps(checks, indent=2) + "\n")
    fields = ("symbol", "ensembl_id", "named_r8_class", "candidate_priority_class", "named_r8_genetics_call", "candidate_priority_genetics_call", "rna_call_unchanged")
    table = "".join('<tr>' + ''.join('<td>' + html.escape(str(row[k])) + '</td>' for k in fields) + '</tr>' for row in changes)
    qualification_link = os.path.relpath(qualification_path, out)
    (out / "index.html").write_text('<meta charset="utf-8"><title>Unadopted Catalog eligibility candidate</title><style>body{font:15px sans-serif;margin:2rem;max-width:1200px}table{border-collapse:collapse}th,td{padding:.4rem;border-bottom:1px solid #ddd;text-align:left}input{padding:.6rem;width:30rem;max-width:90%}</style><h1>Unadopted Catalog eligibility candidate</h1><p><strong>This is a proposed alternate view. Named r8 classifications remain unchanged.</strong> The current priority rule changes 188 genetics evidence states and 180 gene classes. All 28,233 identities and all RNA calls, directions, effects, adjusted P values and missingness are retained. The original class and genetics state remain beside the candidate values.</p><p>Candidate counts: 13,731 untested; 12,958 unresolved; 1,118 established-state associated; 400 genetically anchored; 26 supported in both domains.</p><p>Rule: positive &gt; any untestable study &gt; tested below threshold &gt; no SuSiE result. This proposal does not apply the separate post hoc LD/prior sensitivities and does not qualify unscreened genes as passing those diagnostics.</p><ul><li><a href="catalog_priority_eligibility_candidate.tsv">All 28,233 candidate records, with named r8 and RNA fields preserved</a></li><li><a href="eligibility_changes_only.tsv">188 changed genetics evidence states</a></li><li><a href="' + html.escape(qualification_link) + '">Separate native genetics qualification and sensitivity annotations</a></li><li><a href="checks.json">Independent preservation checks</a></li></ul><input id="q" type="search" placeholder="Filter changed genes by symbol or Ensembl ID" oninput="filterRows()"><table id="changes"><thead><tr>' + ''.join('<th>'+html.escape(k)+'</th>' for k in fields) + '</tr></thead><tbody>' + table + '</tbody></table><script>function filterRows(){const q=document.getElementById("q").value.toLowerCase();document.querySelectorAll("#changes tbody tr").forEach(r=>r.hidden=!r.textContent.toLowerCase().includes(q));}</script>')
    # Update this new review output only; historical Catalog interfaces are untouched.
    index = current / "index.html"
    text = index.read_text()
    links = '<li><a href="' + html.escape(os.path.relpath(out / "index.html", current)) + '">Separate unadopted eligibility-corrected Catalog view</a></li><li><a href="' + html.escape(os.path.relpath(qualification_path, current)) + '">Native genetics qualification and prior/LD annotations</a></li>'
    if "Separate unadopted eligibility-corrected Catalog view" not in text:
        index.write_text(text.replace('</ul>', links + '</ul>', 1))
    print(json.dumps(checks))


if __name__ == "__main__":
    main()
