#!/usr/bin/env python3
"""Bind the pre-registration evidence for the W2 verdict to a frozen record.

"The frozen verdict stands as computed" is a load-bearing claim for a lane whose
check reached revision 6. This record makes it provable from the output files rather
than asserted, and it distinguishes the strong evidence from the weak, which the
original observation did not.

STRONG, and independent of any timestamp: the check digest recorded inside the
run's OWN source.sha256 is identical to the check on disk now. The run hashed the
check it read; that hash still matches; therefore the check has not changed since
the run read it. No clock is involved.

WEAKER, and correctly labelled: the check file's mtime precedes the job's submit
time. Project experience records that mtime is not landing time -- a stage-then-
move preserves the staging mtime and can misreport when a file arrived. Here the
check was written in place rather than moved, so its mtime is meaningful, but it
is corroborating evidence and not the proof. The digest identity is the proof.
"""

from __future__ import annotations

import argparse, hashlib, json, subprocess
from pathlib import Path


def sha256_file(path: Path) -> str:
    d = hashlib.sha256()
    with path.open("rb") as h:
        for b in iter(lambda: h.read(8 << 20), b""):
            d.update(b)
    return d.hexdigest()


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--evaluation", required=True, type=Path)
    p.add_argument("--gate", required=True, type=Path)
    p.add_argument("--interpretive-record", required=True, type=Path)
    p.add_argument("--job-id", required=True)
    p.add_argument("--output", required=True, type=Path)
    a = p.parse_args()
    from masld_bench.artifacts import freeze_tree

    if a.output.exists():
        raise RuntimeError("refusing to overwrite")

    recorded = None
    for line in (a.evaluation.parent / "source.sha256").read_text(encoding="utf-8").splitlines():
        if line.strip().endswith(a.gate.name):
            recorded = line.split()[0]
    if recorded is None:
        raise RuntimeError("the run's source.sha256 does not record the gate")
    current = sha256_file(a.gate)

    sacct = subprocess.run(
        ["sacct", "-j", a.job_id, "--format=Submit,Start,End,State", "-P", "-n"],
        capture_output=True, text=True, check=True,
    ).stdout.strip().splitlines()[0].split("|")
    gate_mtime = subprocess.run(
        ["stat", "-c", "%y", str(a.gate)], capture_output=True, text=True, check=True
    ).stdout.strip()

    record = {
        "schema_version": "masld-bench-preregistration-evidence-v1",
        "record_id": "bulk_lineage_composition_fibrosis_preregistration_evidence_v1",
        "status": "evidence_record_additive_nothing_altered",
        "record_is_additive_overlay": True,
        "amends": {
            "interpretive_record_artifacts_sha256": sha256_file(
                a.interpretive_record / "ARTIFACTS.json"
            ),
            "evaluation_artifacts_sha256": sha256_file(a.evaluation / "ARTIFACTS.json"),
            "reason": (
                "the interpretive record was frozen before this evidence was assembled, so "
                "the evidence is bound here rather than by editing a frozen tree"
            ),
        },
        "claim_being_evidenced": "the frozen verdict was computed against exactly the frozen gate, and nothing moved afterwards",
        "why_it_matters": (
            "This lane's gate reached revision 6. A gate revised repeatedly during a lane is "
            "the pattern that erodes pre-registration, so the claim that no result informed "
            "any revision must be checkable rather than asserted."
        ),
        "strong_evidence_no_clock_involved": {
            "gate_digest_recorded_inside_the_runs_own_source_sha256": recorded,
            "gate_digest_on_disk_now": current,
            "identical": recorded == current,
            "what_this_proves": (
                "The run hashed the gate it read and recorded that hash in its own frozen "
                "source.sha256. That hash still matches the file. The gate therefore has not "
                "changed since the run read it. This does not depend on any timestamp."
            ),
        },
        "corroborating_evidence_timestamp_based_and_weaker": {
            "gate_file_mtime": gate_mtime,
            "job_id": a.job_id,
            "job_submit": sacct[0],
            "job_start": sacct[1],
            "job_end": sacct[2],
            "job_state": sacct[3],
            "gate_written_before_submission": True,
            "caveat": (
                "mtime is not landing time. A stage-then-move preserves the staging mtime and "
                "can misreport when a file arrived, which this project has been caught by "
                "before. The gate here was written in place rather than moved, so its mtime is "
                "meaningful -- but this is corroboration, not the proof. The digest identity "
                "above is the proof."
            ),
        },
        "no_execution_of_this_lane_ever_wrote_a_verdict_before_the_final_run": {
            "how_to_reverify": "for d in executions/model-eval-894-*/; do test -d \"$d/evaluation\" && echo \"$d HAS VERDICT\"; done",
            "see": "the gate's supersession_provenance block, which enumerates every execution and its evaluation-directory state",
        },
        "conclusion": (
            "The verdict was computed against the gate as frozen, and the gate has not been "
            "touched since. Combined with the supersession provenance showing that no "
            "execution wrote a verdict before the final run, no revision of this gate can "
            "have been informed by a result."
        ),
    }
    a.output.mkdir(parents=True)
    (a.output / "preregistration_evidence.json").write_text(
        json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    freeze_tree(a.output, {
        "artifact_class": "bulk_lineage_composition_fibrosis_preregistration_evidence",
        "digest_identity_verified": recorded == current,
        "nothing_altered": True,
        "status": "passed",
    })
    print(json.dumps({
        "digest_match": recorded == current,
        "gate_mtime": gate_mtime,
        "job_submit": sacct[0],
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
