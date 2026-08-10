#!/usr/bin/env python3
"""Freeze evidence-class producer code and covariate sources before execution."""

from __future__ import annotations

import datetime as dt
import json
from pathlib import Path

from public_functional_common import CANDIDATE_ROOT, PROJECT_ROOT, SCRIPT_ROOT, atomic_write_text, require_sealed, sha256_file, stable_json_sha256, write_tsv


GTF = Path("/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz")
GTF_SHA256 = "73bbbbd6eb2f114d1536f7cbf2339a653e1edc889b0cbe78992e92560b5aaba9"


def main() -> None:
    seal = require_sealed()
    destination = CANDIDATE_ROOT / "EVIDENCE_CLASS_CODE_READY.json"
    if destination.exists():
        raise RuntimeError(f"Evidence-class execution already frozen: {destination}")
    if sha256_file(GTF) != GTF_SHA256:
        raise RuntimeError("GENCODE v49 GTF drift")
    paths = [SCRIPT_ROOT / "functional_analysis_common.R", SCRIPT_ROOT / "07_run_evidence_classes.R", SCRIPT_ROOT / "run_evidence_classes.sbatch"]
    rows = [{"path": str(path.relative_to(PROJECT_ROOT)), "sha256": sha256_file(path), "size_bytes": path.stat().st_size} for path in paths]
    rows.append({"path": str(GTF), "sha256": GTF_SHA256, "size_bytes": GTF.stat().st_size})
    contract = {
        "freeze_id": "public-functional-evidence-class-v1",
        "frozen_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "parent_specification_sha256": seal["specification_sha256"],
        "files": rows,
        "primary_method": "sample-level limma camera target class versus neither, with empirical per-set inter-gene correlation",
        "alignment": "induction=sign(canonical_bulk_logFC)*effect; reversal=-sign(canonical_bulk_logFC)*effect",
        "sensitivity": "deterministic no-replacement nearest matching on assay mean, assay variance, and GENCODE-v49 gene length; 10,000 paired sign flips",
        "multiplicity": "BH over disease_state_only, genetic_only, and convergent separately per dataset contrast and method",
    }
    contract["code_bundle_sha256"] = stable_json_sha256(contract)
    write_tsv(CANDIDATE_ROOT / "evidence_class_code_manifest.tsv", rows, ["path", "sha256", "size_bytes"])
    atomic_write_text(destination, json.dumps(contract, indent=2, sort_keys=True) + "\n")
    print(json.dumps(contract, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
