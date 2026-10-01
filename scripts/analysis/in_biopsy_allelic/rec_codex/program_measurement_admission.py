"""Metadata-only admission of a filtered RNA fixture to fixed program targets."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import sys

ROOT = Path(__file__).resolve().parents[4]
MEMBERS = ROOT / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_membership_v2.tsv"
FIXTURE = ROOT / "Analysis/MASLD_Model_Benchmark/executions/model-data-079-21083342/model_input"
PANEL = ROOT / "RNA-seq/results/cross_assay_modules/cam-v2-20260912T130258Z/universe/measured_cosmx_govaere.txt"
EXPECTED = {
    "genes": "8d3fa19fa6495a213d6b10ab350ed0cae1c2bf64edc36c758ea9421c6a13972e",
    "receipt": "cd20f55f886fd3a3ca9b8a9584c83d5815b32b6a86a6206eb0af99a92a698cc3",
    "membership": "feec3fc9ccaaffaa9abc1ac5d593cc06845ed8140460c8873a51bc8652d7336b",
    "panel": "2d973d7e0ba08a2c3aaf92d2c7fb296b1d4924eb294be9cc9e420a48af1e4311",
    "participants": "01f63880e78f179a5a8670434fc3d33fd92f615464dc4efd8d80ee2598dac52c",
}
SPLIT_HASHES = {
    "fit": "07a7ee166ed8aaa79dcc63cf1552d7b46c897b8649601f145f36b678786a3b7e",
    "calibration": "b157d378921c8154d0639993f49020d7f2a3e50031f31cabe321d36f9146d830",
    "internal_test": "9efe3fc79f2873ddd8e5ff68e52dab99f374ba83598475586315a3122d008448",
}
CONFIRMED = "gencode_v49_unique_symbol_confirmed"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields):
    with path.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main(out):
    require(os.environ.get("SLURM_JOB_ID"), "Run metadata admission in a SLURM allocation.")
    paths = {"membership": MEMBERS, "panel": PANEL,
             "genes": FIXTURE / "candidate_gene_axis.tsv",
             "participants": FIXTURE / "participant_axis.tsv", "receipt": FIXTURE / "receipt.json"}
    hashes = {name: sha(path) for name, path in paths.items()}
    for name, expected in EXPECTED.items():
        require(hashes[name] == expected, f"Frozen {name} hash changed; do not silently redefine admission.")
    members, genes, people = read_tsv(MEMBERS), read_tsv(paths["genes"]), read_tsv(paths["participants"])
    receipt = json.loads(paths["receipt"].read_text())
    panel = PANEL.read_text().splitlines()
    require(len(panel) == len(set(panel)) == 968 and all(panel), "Panel must have 968 unique names.")
    require(len(genes) == 14078, "Unexpected filtered fixture feature count.")
    require([int(g["feature_index"]) for g in genes] == list(range(len(genes))), "Fixture feature index is not contiguous.")
    require(len({g["stable_gene_id"] for g in genes}) == len(genes), "Repeated stable-gene columns.")
    require(all(g["gencode_v49_gene_name"] for g in genes), "Fixture has unnamed genes.")
    require(len(people) == 109 and len({p["row_id"] for p in people}) == 109, "109 unique participant rows required.")
    require([int(p["row_index"]) for p in people] == list(range(109)), "Participant index is not contiguous.")
    require(all(p["cohort_family_id"] == "gse268273_imid_masld" and
                p["rna_observation_state"] == "observed" and
                int(p["raw_technical_run_count"]) > 0 and
                p["quantification_measurement"] == "RSEM_expected_count_raw_count_scale" for p in people),
            "Participant measurement/unit metadata changed.")
    require(receipt["shape"] == [109, 14078] and receipt["participants"] == 109 and
            receipt["measurement"] == "RSEM_expected_count_raw_count_scale" and
            receipt["normalization_applied"] is False and receipt["rna_observed_mask_all_true"] is True,
            "Fixture receipt does not describe the admitted raw-count metadata.")
    by_symbol = defaultdict(list)
    for gene in genes:
        by_symbol[gene["gencode_v49_gene_name"]].append(gene)
    groups = defaultdict(list)
    for member in members:
        groups[member["program_uid"]].append(member)
    require(len(groups) == 117 and len(members) == 7093, "Frozen 117-program membership shape changed.")
    coverage, member_rows, admitted_weights = [], [], []
    for uid, group in sorted(groups.items()):
        require(len({m["membership_sha256"] for m in group}) == 1, f"Inconsistent membership hash: {uid}")
        require(len({m["source_gene"] for m in group}) == len(group), f"Repeated source member: {uid}")
        weights = [float(m["original_l1_weight"]) for m in group]
        raw = [abs(float(m["canonical_weight_text"])) for m in group]
        require(all(math.isfinite(w) and w >= 0 for w in weights) and
                all(math.isfinite(w) for w in raw) and math.fsum(raw) > 0, f"Invalid source weights: {uid}")
        require(math.isclose(math.fsum(weights), 1, abs_tol=1e-10, rel_tol=0), f"Source L1 weights do not sum to one: {uid}")
        require(all(math.isclose(w, r / math.fsum(raw), abs_tol=1e-12, rel_tol=1e-9)
                    for w, r in zip(weights, raw)), f"Original weights differ from normalized source weights: {uid}")
        counts = Counter()
        program_rows = []
        for member in group:
            symbol = member["mapped_symbol"]
            unconfirmed = member["mapped_symbol_status"] != CONFIRMED
            absent = symbol not in by_symbol
            duplicate = len(by_symbol.get(symbol, [])) > 1
            reasons = [name for name, failed in [("mapping_unconfirmed", unconfirmed),
                       ("absent_filtered_fixture", absent), ("duplicated_count_symbol", duplicate)] if failed]
            counts.update(reasons)
            measured = symbol in panel and not reasons
            row = {**member, "mapping_unconfirmed": int(unconfirmed), "absent_filtered_fixture": int(absent),
                   "duplicated_count_symbol": int(duplicate), "fixture_column_count": len(by_symbol.get(symbol, [])),
                   "fixture_feature_indices": ",".join(g["feature_index"] for g in by_symbol.get(symbol, [])),
                   "admitted_panel_member": int(measured), "exclusion_reasons": ";".join(reasons)}
            program_rows.append(row)
        complete = not counts
        for row in program_rows:
            row["program_complete_target"] = int(complete)
        member_rows.extend(program_rows)
        if complete:
            admitted_weights.extend(program_rows)
        coverage.append({"program_uid": uid, "source_members": len(group), "source_l1_sum": math.fsum(weights),
                         "complete_fixture_target": int(complete), "mapping_unconfirmed_members": counts["mapping_unconfirmed"],
                         "absent_filtered_fixture_members": counts["absent_filtered_fixture"],
                         "duplicated_count_symbol_members": counts["duplicated_count_symbol"],
                         "excluded_members_union": sum(bool(r["exclusion_reasons"]) for r in program_rows),
                         "admitted_panel_members": sum(r["admitted_panel_member"] for r in program_rows),
                         "retained_original_l1_weight": math.fsum(float(r["original_l1_weight"]) for r in program_rows if r["admitted_panel_member"]),
                         "target_state": "complete_filtered_fixture" if complete else "incomplete_filtered_fixture"})
    panel_rows = []
    for symbol in panel:
        columns = by_symbol.get(symbol, [])
        panel_rows.append({"symbol": symbol, "fixture_column_count": len(columns),
                           "fixture_feature_indices": ",".join(g["feature_index"] for g in columns),
                           "admitted": int(len(columns) == 1),
                           "exclusion_reason": "" if len(columns) == 1 else
                           "absent_filtered_fixture" if not columns else "duplicated_count_symbol"})
    ordered = sorted(people, key=lambda p: (hashlib.sha256(("20260930\tmissing117\t" + p["row_id"]).encode()).hexdigest(), p["row_id"]))
    allocation = []
    for i, person in enumerate(ordered):
        role = "fit" if i < 69 else "calibration" if i < 89 else "internal_test"
        allocation.append({**person, "split_order": i, "role": role, "inner_fold": i % 5 if role == "fit" else ""})
    split_hashes = {}
    for role, expected in SPLIT_HASHES.items():
        ids = [p["row_id"] for p in allocation if p["role"] == role]
        split_hashes[role] = hashlib.sha256(("\n".join(ids) + "\n").encode()).hexdigest()
        require(split_hashes[role] == expected, f"Fixed {role} allocation changed.")
    require(Counter(p["inner_fold"] for p in allocation if p["role"] == "fit") == {0: 14, 1: 14, 2: 14, 3: 14, 4: 13}, "Inner fold counts changed.")
    out.mkdir(parents=True, exist_ok=False)
    for name, rows in [("program_coverage.tsv", coverage), ("member_admission.tsv", member_rows),
                       ("admitted_original_weights.tsv", admitted_weights), ("panel_axis_admission.tsv", panel_rows),
                       ("participant_allocation.tsv", allocation)]:
        write_tsv(out / name, rows, list(rows[0]) if rows else list(member_rows[0]))
    summary = {"schema_version": "fixture-program-admission-v1", "scope": "filtered_fixture_metadata_only",
               "source_files": {name: {"path": str(path), "sha256": hashes[name]} for name, path in paths.items()},
               "script_sha256": sha(__file__), "python": sys.version, "platform": platform.platform(),
               "job_id": os.environ["SLURM_JOB_ID"], "programs": 117,
               "complete_targets": sum(p["complete_fixture_target"] for p in coverage),
               "source_columns": len(genes), "distinct_symbols": len(by_symbol),
               "duplicate_symbols": sorted(k for k, v in by_symbol.items() if len(v) > 1),
               "panel_names": len(panel), "unambiguous_panel_inputs": sum(p["admitted"] for p in panel_rows),
               "reason_counts_are_overlapping": True,
               "member_reason_counts": {key: sum(p[key] for p in coverage) for key in
                   ["mapping_unconfirmed_members", "absent_filtered_fixture_members", "duplicated_count_symbol_members"]},
               "split_counts": dict(Counter(p["role"] for p in allocation)), "ordered_split_sha256": split_hashes,
               "inferential_unit": "participant metadata row; universal cross-study person disjointness unverified",
               "technical_runs_kept_together": True, "weights_renormalized": False,
               "expression_arrays_read": False, "outcomes_read": False, "fit_performed": False,
               "source_assay_absence_inferred": False,
               "limits": "Fixture gaps do not establish biological absence; duplicate symbols refused, not summed; full-source axis assessment separate."}
    summary["output_sha256"] = {path.name: sha(path) for path in sorted(out.glob("*.tsv"))}
    (out / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"complete_fixture_targets": summary["complete_targets"],
                      "unambiguous_panel_inputs": summary["unambiguous_panel_inputs"], "output": str(out)}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    main(parser.parse_args().output)
