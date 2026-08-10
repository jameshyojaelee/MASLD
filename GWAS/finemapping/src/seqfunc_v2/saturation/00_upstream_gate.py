#!/usr/bin/env python3
"""Hard gate uniform35 completeness and adult-hepatocyte model validity.

The gate deliberately accepts explicit contract paths because the ChromBPNet
branch is produced independently. It never infers PASS from file presence
alone and emits a scorer-ready one-row-per-fold model manifest.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

from common import CHROMBPNET, FOLD_GROUPS, OUT, UNIFORM, atomic_json, read_json, sha256, truthy


def pick(obj, paths):
    for path in paths:
        cur = obj
        try:
            for key in path.split("."):
                cur = cur[key]
        except (KeyError, TypeError):
            continue
        return cur
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--uniform-audit", type=Path, default=UNIFORM / "audit/audit_verdict.json")
    ap.add_argument("--chrombpnet-gate", type=Path, default=CHROMBPNET / "gate_verdict.json")
    ap.add_argument("--model-manifest", type=Path, default=CHROMBPNET / "fold_manifest.tsv")
    ap.add_argument("--run-contract", type=Path, default=CHROMBPNET / "run_contract.json")
    ap.add_argument("--out-dir", type=Path, default=OUT / "gates")
    args = ap.parse_args()
    checks = []

    def check(name, passed, detail):
        checks.append({"check": name, "pass": bool(passed), "detail": str(detail)})

    if not args.uniform_audit.is_file():
        check("uniform35_audit_exists", False, args.uniform_audit)
    else:
        u = read_json(args.uniform_audit)
        check("uniform35_pipeline_complete", truthy(u.get("pipeline_complete")), u.get("pipeline_complete"))
        # Global release coverage is expected to remain false when PolyFun
        # blocks are absent/truncated. Saturation is restricted downstream to
        # primary_eligible reliable loci, so this is recorded, not overridden.
        check("uniform35_release_gate_recorded_not_required", "release_gate_pass" in u,
              f"release_gate_pass={u.get('release_gate_pass')}")
        check("uniform35_canonical_firewall", u.get("canonical_outputs_mutated") is False,
              u.get("canonical_outputs_mutated"))
        locus_manifest = UNIFORM / "config/locus_manifest.tsv"
        locus_summary = UNIFORM / "aggregate/locus_summary.tsv"
        if not locus_manifest.is_file() or not locus_summary.is_file():
            check("uniform35_terminal_accounting_files", False, f"{locus_manifest};{locus_summary}")
        else:
            with locus_manifest.open() as h: expected = list(csv.DictReader(h, delimiter="\t"))
            with locus_summary.open() as h: observed = list(csv.DictReader(h, delimiter="\t"))
            terminal = {"completed","nonconverged","insufficient_ld","allele_failure","ld_failure","model_failure"}
            exact = (len(expected) == len(observed) and
                     {x["locus_index"] for x in expected} == {x["locus_index"] for x in observed} and
                     all(x.get("status") in terminal for x in observed))
            check("uniform35_exact_terminal_accounting", exact,
                  f"expected={len(expected)} observed={len(observed)} terminal={sum(x.get('status') in terminal for x in observed)}")

            # Terminal is not the same as successful.  The accounting check above
            # passes when every locus reached SOME terminal status, including
            # mass failure -- which is how a run that modelled 17 of 1,950 loci
            # unlocked GPU scoring.  Gate on yield explicitly, as defence in
            # depth behind the audit's own gates (04_audit.R).
            MIN_COMPLETION_RATE = 0.90
            modelled = sum(x.get("status") in {"completed", "nonconverged"} for x in observed)
            rate = (modelled / len(observed)) if observed else 0.0
            check("uniform35_yield", rate >= MIN_COMPLETION_RATE,
                  f"{modelled}/{len(observed)} = {rate:.3f} (min {MIN_COMPLETION_RATE})")
            tier1 = [x for x in observed if str(x.get("tier")) == "1"]
            tier1_done = sum(x.get("status") == "completed" for x in tier1)
            check("uniform35_tier1_yield", len(tier1) > 0 and tier1_done >= 1,
                  f"{tier1_done} completed of {len(tier1)} tier-1 loci")

    if not args.chrombpnet_gate.is_file():
        check("chrombpnet_gate_exists", False, args.chrombpnet_gate)
    else:
        c = read_json(args.chrombpnet_gate)
        internal = pick(c, ["fold0_internal_qc_pass", "internal_gate_pass", "internal.pass"])
        currin = pick(c, ["currin_external_calibration_pass", "currin_external_gate_pass", "external_currin.pass"])
        unlocked = c.get("saturation_unlocked")
        check("chrombpnet_internal_gate", truthy(internal), internal)
        check("currin_external_gate", truthy(currin), currin)
        check("producer_saturation_unlocked", truthy(unlocked) and str(c.get("status","")).lower() == "pass",
              f"unlocked={unlocked};status={c.get('status')}")

    reference_fasta = None
    chrom_sizes = None
    if not args.run_contract.is_file():
        check("chrombpnet_run_contract_exists", False, args.run_contract)
    else:
        rc = read_json(args.run_contract)
        check("calibration_training_prohibited", rc.get("calibration_labels_used_for_training") is False,
              rc.get("calibration_labels_used_for_training"))
        check("donor_heldout_claim_prohibited", rc.get("donor_heldout_claim_supported") is False,
              rc.get("donor_heldout_claim_supported"))
        check("adult_hepatocyte_scope", rc.get("model_cell_type") == "adult_liver_hepatocyte",
              rc.get("model_cell_type"))
        reference_fasta = Path(rc.get("reference_fasta", ""))
        chrom_sizes = CHROMBPNET / "GRCh38.autosomes.chrom.sizes"
        check("producer_reference_fasta_exists", reference_fasta.is_file() and reference_fasta.stat().st_size > 0,
              reference_fasta)
        check("producer_chrom_sizes_exists", chrom_sizes.is_file() and chrom_sizes.stat().st_size > 0,
              chrom_sizes)

    models = []
    if not args.model_manifest.is_file():
        check("model_manifest_exists", False, args.model_manifest)
    else:
        with args.model_manifest.open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        check("five_model_rows", len(rows) == 5, len(rows))
        by_fold = {("fold" + row["fold"] if row.get("fold","").isdigit() else row.get("fold_id")): row for row in rows}
        check("exact_fold_ids", set(by_fold) == set(FOLD_GROUPS), sorted(x for x in by_fold if x))
        for fold, chroms in FOLD_GROUPS.items():
            row = by_fold.get(fold, {})
            observed = sorted(filter(None, row.get("test_chromosomes", "").replace(",",";").split(";")))
            check(f"{fold}_exact_test_chromosomes", observed == sorted(chroms), observed)
            model = Path(row.get("model_path") or CHROMBPNET/fold/"chrombpnet_full/models/chrombpnet_nobias.h5")
            peaks = Path(row.get("peaks_path") or CHROMBPNET/"hepatocyte.idr_pooled_summit.narrowPeak")
            fold_json = Path(row.get("fold_json") or CHROMBPNET/fold/f"{fold}.json")
            for label, path in (("model", model), ("peaks", peaks), ("fold_json", fold_json)):
                check(f"{fold}_{label}_exists", path.is_file() and path.stat().st_size > 0, path)
            reference_ready = (reference_fasta is not None and chrom_sizes is not None and
                               reference_fasta.is_file() and chrom_sizes.is_file())
            if reference_ready and all(p.is_file() and p.stat().st_size > 0 for p in (model, peaks, fold_json)):
                fobj = read_json(fold_json)
                test = sorted(fobj.get("test", fobj.get("test_chromosomes", [])))
                check(f"{fold}_json_matches", test == sorted(chroms), test)
                iqc_path = CHROMBPNET / fold / "internal_qc.json"
                if not iqc_path.is_file():
                    check(f"{fold}_internal_qc_exists", False, iqc_path)
                else:
                    iqc = read_json(iqc_path)
                    check(f"{fold}_internal_qc_pass", truthy(iqc.get("internal_qc_pass")), iqc.get("internal_qc_pass"))
                models.append({
                    "fold_id": fold, "test_chromosomes": ";".join(chroms),
                    "model_path": str(model.resolve()), "model_sha256": sha256(model),
                    "peaks_path": str(peaks.resolve()), "peaks_sha256": sha256(peaks),
                    "fold_json": str(fold_json.resolve()), "fold_json_sha256": sha256(fold_json),
                    "reference_fasta": str(reference_fasta.resolve()),
                    "chrom_sizes_path": str(chrom_sizes.resolve()), "chrom_sizes_sha256": sha256(chrom_sizes),
                })

    passed = bool(checks) and all(x["pass"] for x in checks)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    with (args.out_dir / "upstream_checks.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["check", "pass", "detail"], delimiter="\t")
        writer.writeheader(); writer.writerows(checks)
    if passed:
        with (args.out_dir / "model_manifest.locked.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(models[0]), delimiter="\t")
            writer.writeheader(); writer.writerows(sorted(models, key=lambda x: x["fold_id"]))
    atomic_json({
        "status": "PASS" if passed else "BLOCKED", "saturation_submission_allowed": passed,
        "uniform_audit": str(args.uniform_audit), "chrombpnet_gate": str(args.chrombpnet_gate),
        "source_model_manifest": str(args.model_manifest), "chrombpnet_run_contract": str(args.run_contract), "checks": checks,
        "apply_only_firewall": True, "canonical_outputs_mutated": False,
    }, args.out_dir / "upstream_gate.json")
    if not passed:
        raise SystemExit("Saturation hard gate BLOCKED; see upstream_checks.tsv")
    print("Saturation upstream gate PASS")


if __name__ == "__main__":
    main()
