#!/usr/bin/env python3
"""Assemble the apply-only SeqFunc v2 release gate without promoting outputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(os.environ.get("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"))
# No default release root.  This used to default to 2026-07-13-r1, so running
# the gate bare rewrote that FROZEN release's verdict in place -- and it is the
# older, FAILed release, which docs/STATUS.md cites as a re-entry
# check.  --release-root is now required.
RELEASE = None


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path):
    with path.open() as handle:
        return json.load(handle)


def atomic_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f".tmp.{os.getpid()}")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    os.replace(tmp, path)


def atomic_tsv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f".tmp.{os.getpid()}")
    fields = list(rows[0]) if rows else ["stage", "state", "detail"]
    with tmp.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp, path)


def boolean_from_gate(value) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.upper() == "PASS"
    return False


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--release-root", type=Path, default=RELEASE, required=True)
    parser.add_argument("--force", action="store_true",
                        help="re-run the gate over a release that already has a verdict")
    args = parser.parse_args()
    release = args.release_root.resolve()
    # canonical_firewall_after.tsv, release_checks.tsv and release_verdict.json
    # are rewritten IN PLACE by this script.  Re-running it over a release that
    # already has a verdict destroys that frozen record (and restamps its
    # timestamp), so refuse unless explicitly forced.
    existing_verdict = release / "release_verdict.json"
    if existing_verdict.exists() and not args.force:
        raise SystemExit(
            f"{existing_verdict} already exists -- re-running would overwrite a frozen verdict. "
            "Pass --force only if you intend to replace it."
        )
    contract = read_json(release / "run_contract.json")
    checks: list[dict] = []

    def add(stage: str, state: str, detail: str, path: Path | None = None) -> None:
        if state not in {"PASS", "FAIL", "PENDING", "REPORTED"}:
            raise ValueError(state)
        checks.append({
            "stage": stage,
            "state": state,
            "detail": detail,
            "path": str(path.resolve()) if path is not None else "",
        })

    # P0 frozen contracts.
    studies = list(csv.DictReader((release / "study_manifest.tsv").open(), delimiter="\t"))
    folds = list(csv.DictReader((release / "fold_manifest.tsv").open(), delimiter="\t"))
    add("frozen_contracts", "PASS" if len(studies) == 35 and len(folds) == 5 else "FAIL",
        f"studies={len(studies)} folds={len(folds)}", release / "run_contract.json")

    # Recompute every protected hash. This is always a hard gate.
    before = list(csv.DictReader((release / "canonical_firewall_before.tsv").open(), delimiter="\t"))
    after = []
    firewall_ok = True
    for row in before:
        path = Path(row["path"])
        observed = sha256(path) if path.is_file() else "MISSING"
        same = observed == row["sha256"]
        firewall_ok &= same
        after.append({"path": str(path), "before_sha256": row["sha256"], "after_sha256": observed,
                      "unchanged": str(same).upper()})
    atomic_tsv(release / "canonical_firewall_after.tsv", after)
    add("canonical_firewall", "PASS" if firewall_ok else "FAIL",
        f"unchanged={sum(r['unchanged'] == 'TRUE' for r in after)}/{len(after)}",
        release / "canonical_firewall_after.tsv")

    currin = ROOT / "GWAS/finemapping/results/seqfunc/chrombpnet_caqtl/v2_truth/truth_contract.json"
    if currin.is_file():
        c = read_json(currin)
        ok = (c.get("n_positive_1kb") == 35361 and c.get("n_source_background") == 125689
              and c.get("one_mb_signed_use") is False)
        add("currin_truth", "PASS" if ok else "FAIL",
            f"primary={c.get('n_positive_1kb')} background={c.get('n_source_background')}", currin)
    else:
        add("currin_truth", "PENDING", "truth contract absent", currin)

    finemap = Path(contract["paths"]["finemap"]) / "audit/audit_verdict.json"
    if finemap.is_file():
        f = read_json(finemap)
        complete = f.get("pipeline_complete") is True
        passed = f.get("release_gate_pass") is True
        mutated = f.get("canonical_outputs_mutated") is not False   # True or missing => integrity break
        hard = list(f.get("hard_failures") or [])
        # A pipeline that ran to completion with no hard failure and no canonical
        # mutation, but whose release gate declines to promote for a documented
        # coverage reason (e.g. all_gws_represented under strict PolyFun LD), is a
        # REPORTED honest negative -- not a release FAIL. Integrity breaks still FAIL.
        if mutated or hard:
            state = "FAIL"
        elif complete and passed:
            state = "PASS"
        elif complete:
            state = "REPORTED"
        else:
            state = "PENDING"
        add("uniform35", state,
            f"pipeline_complete={complete} release_gate_pass={passed} "
            f"release_failures={f.get('release_failures')} hard_failures={len(hard)}", finemap)
    else:
        add("uniform35", "PENDING", "fine-mapping audit absent", finemap)

    splicing_dir = Path(contract["paths"]["splicing"])
    splicing = splicing_dir / "RELEASE_GATE.tsv"
    perm_gate = splicing_dir / "permutations/permutation_gate.json"
    if splicing.is_file():
        ok = splicing.read_text().split("\t", 1)[0].strip() == "PASS"
        add("leafcutter", "PASS" if ok else "FAIL", splicing.read_text().strip(), splicing)
    elif perm_gate.is_file() and read_json(perm_gate).get("final_fail") is True:
        # LeafCutter clustering completed but the permutation null failed calibration
        # (documented honest negative); the release gate was honestly not emitted.
        p = read_json(perm_gate)
        add("leafcutter", "REPORTED",
            f"honest permutation-null failure final_fail=True "
            f"median_null_lambda={p.get('median_null_lambda')} n_permutations={p.get('n_permutations')}",
            perm_gate)
    else:
        add("leafcutter", "PENDING", "release gate absent", splicing)

    mpra = Path(contract["paths"]["mpra"]) / "gate_verdict.json"
    if mpra.is_file():
        m = read_json(mpra)
        source_ok = boolean_from_gate(m.get("source_reproduction_pass"))
        source_done = boolean_from_gate(m.get("source_reproduction_complete"))
        sensitivity_done = boolean_from_gate(
            m.get("mpranalyze_complete", m.get("mpranalyze_completion"))
        )
        mpra_state = "PASS" if source_ok and sensitivity_done else (
            "FAIL" if source_done and not source_ok else "PENDING"
        )
        add("hu_mpra", mpra_state,
            f"source_reproduction={source_ok} mpranalyze_complete={sensitivity_done}", mpra)
    else:
        add("hu_mpra", "PENDING", "branch gate absent", mpra)

    chrom = Path(contract["paths"]["chrombpnet"]) / "gate_verdict.json"
    saturation_unlocked = False
    if chrom.is_file():
        c = read_json(chrom)
        internal = boolean_from_gate(c.get("fold0_internal_qc_pass"))
        external = boolean_from_gate(c.get("currin_external_calibration_pass"))
        saturation_unlocked = boolean_from_gate(c.get("saturation_unlocked")) and internal and external
        awaiting = str(c.get("status", "")).startswith("awaiting")
        add("chrombpnet", "PASS" if saturation_unlocked else ("PENDING" if awaiting else "FAIL"),
            f"fold0_internal={internal} currin_external={external} saturation_unlocked={saturation_unlocked}", chrom)
    else:
        add("chrombpnet", "PENDING", "branch gate absent", chrom)

    saturation = Path(contract["paths"]["saturation"]) / "gate_verdict.json"
    if saturation.is_file():
        s = read_json(saturation)
        ok = saturation_unlocked and boolean_from_gate(s.get("complete"))
        add("saturation", "PASS" if ok else "FAIL",
            f"complete={s.get('complete')} upstream_unlocked={saturation_unlocked}", saturation)
    else:
        add("saturation", "PENDING", "gated downstream branch not yet complete", saturation)

    states = [row["state"] for row in checks]
    if "FAIL" in states:
        overall = "FAIL"
    elif "PENDING" in states:
        overall = "PENDING"
    elif all(x == "PASS" for x in states):
        overall = "PASS"
    else:
        # only PASS + REPORTED remain: every branch reached a defensible terminal
        # state (pass or a documented honest negative), so the apply-only release
        # is complete without claiming every sensitivity branch passed.
        overall = "COMPLETE_WITH_HONEST_NEGATIVES"
    verdict = {
        "release_id": release.name,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "overall": overall,
        "apply_only": True,
        "canonical_promotion": False,
        "saturation_unlocked": saturation_unlocked,
        "checks": checks,
    }
    atomic_tsv(release / "release_checks.tsv", checks)
    atomic_json(release / "release_verdict.json", verdict)
    print(json.dumps(verdict, indent=2))
    if overall == "FAIL":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
